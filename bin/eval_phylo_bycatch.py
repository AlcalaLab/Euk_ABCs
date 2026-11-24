#!/usr/bin/env python3

import glob, json, math, statistics, subprocess, time

import pandas as pd

from datetime import timedelta
from ete4 import Tree
from pathlib import Path

from bin import ncbi_taxonomy as nt


def check_file_exists(file: str) -> bool:
    return (Path(file).is_file() and Path(file).stat().st_size > 0)


def convert_jplace_for_ete4(tree_string):
        # replaces the {edge_num} tags to NHX format to make compatible with ete4
    fixed_tree = tree_string.replace("{", "[&&NHX:edge_num=").replace("}", "]")

    tree = Tree(fixed_tree) #, parser = 1)

    edge_to_node = {}
    for node in tree.traverse():
        # just care about the nodes that our query fits into
        if "edge_num" in node.props:
            edge_to_node[int(node.props["edge_num"])] = node

    return tree, edge_to_node


def check_json_fields(json_fields: list):
    json_field_dict = {'edge_num': json_fields.index('edge_num'),
        'like_weight_ratio': json_fields.index('like_weight_ratio'),
        'distal_length': json_fields.index('distal_length'),
        'pendant_length': json_fields.index('pendant_length')}

    return json_field_dict


def parse_jplace_file(jplace_file):
    """Loads a .jplace file and returns the ete4 Tree and a dict of query placements."""

    jplace_data = json.load(open(jplace_file))

    tree, edge_to_node = convert_jplace_for_ete4(jplace_data['tree'])

    fields_dict = check_json_fields(jplace_data['fields'])

    by_query = {}

    for place in jplace_data['placements']:
        parsed_list = []
        for p in place['p']:
                parsed_list.append(
                    {'edge_num': int(p[fields_dict['edge_num']]),
                        'lwr': float(p[fields_dict['like_weight_ratio']]),
                        'distal': float(p[fields_dict['distal_length']]),
                        'pendant': float(p[fields_dict['pendant_length']]),
                        'node': edge_to_node[int(p[fields_dict['edge_num']])]})

        names = place.get("n")
        for qname in names:
            by_query[qname] = parsed_list

    return tree, by_query


def assess_gappa_edpl(
        outdir: str,
        jplace_file: str): # -> dict[str, float] | None:
    """replaces the calculate_edpl"""

    tmp_gappa_dir = f'{outdir}Phylogeny_Based/Temp_Gappa/'

    Path(tmp_gappa_dir).mkdir(exist_ok = True, parents = True)

    edpl_csv = f'{jplace_file.rpartition("/")[-1].rstrip('jplace')}'

    edpl_cmd = 'gappa examine edpl ' \
        '--allow-file-overwriting ' \
        f'--jplace-path {jplace_file} ' \
        f'--out-dir {tmp_gappa_dir} ' \
        f'--file-prefix {edpl_csv}'

    if not check_file_exists(f'{edpl_csv}edpl_list.csv'):
        # print(edpl_cmd)
        subprocess.run(
            edpl_cmd,
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

    return f'{tmp_gappa_dir}{edpl_csv}edpl_list.csv'


def load_gappa_edpl(outdir: str,
        jplace_file: str):
    """Parses gappa's edpl.csv into a {query_name: float_edpl} dictionary."""

    edpl_df = pd.read_csv(assess_gappa_edpl(outdir, jplace_file))

    edpl_dict = edpl_df.set_index('Pquery')['EDPL'].to_dict()

    return edpl_dict



def calculate_normalized_entropy(query_placements):
    # calculates the normalized Shannon entropy. This is essentially just placement
    # "uncertainty" (0: certain, 1: no idea)
    if len(query_placements) <= 1:
        return 0.0

    total_lwr = sum(qp['lwr'] for qp in query_placements)
    probs = [qp['lwr'] / total_lwr for qp in query_placements if qp['lwr'] > 0]

    if len(probs) <= 1:
        return 0.0

    raw_entropy = -sum(prob * math.log(prob) for prob in probs)
    max_entropy = math.log(len(probs))

    return raw_entropy / max_entropy



def find_cumulative_lca(tree, query_placements, lwr_cutoff: float = 0.80):
    """Finds the LCA node covering candidate edges up to the cumulative LWR cutoff."""

    sorted_qp = sorted(query_placements, key = lambda x: x['lwr'], reverse = True)

    agg_lwr = 0.0
    target_nodes = []

    for qp in sorted_qp:
        agg_lwr += qp['lwr']

        target_nodes.append(qp['node'])

        if agg_lwr >= lwr_cutoff:
            break

    unique_nodes = list(dict.fromkeys(target_nodes))

    if len(unique_nodes) == 1:
        lca_node = unique_nodes[0]

    elif tree in unique_nodes:
        lca_node = tree

    else:
        try:
            unique_leaves = [leaf.name for target_node in target_nodes for leaf in target_node.leaves()]
            lca_node = tree.common_ancestor(list(set(unique_leaves)))

        except:
            lca_node = tree

    return lca_node, agg_lwr


def query_uncertainty_summary(
        query_name,
        tree,
        query_placements,
        edpl_dict,
        lwr_cutoff: float = 0.80):

    """Combines gappa EDPL with Python entropy and cumulative LCA checks."""
    if not query_placements:
        return None

    sorted_qp = sorted(query_placements, key = lambda x: x['lwr'], reverse = True)

    lca_node, agg_lwr = find_cumulative_lca(tree, sorted_qp, lwr_cutoff)

    uncertainty_summary = {'Top_Edge': sorted_qp[0]['edge_num'],
            'Top_LWR': sorted_qp[0]['lwr'],
            'Total_Placements': len(sorted_qp),
            'EDPL': edpl_dict.get(query_name, 0.0),
            'Normalized_Entropy': calculate_normalized_entropy(sorted_qp),
            'LCA_Node': lca_node,
            'LCA_Agg_LWR': agg_lwr,
            'Mean_Pendant': sum(qp['pendant'] * qp['lwr'] for qp in sorted_qp)}

    return uncertainty_summary


def tree_thresholds(tree, blen_mode: str = 'median'):
    """Derives tree-scaled EDPL and pendant length cutoffs from reference branch lengths."""
    branch_lengths = [node.dist for node in tree.traverse() if not node.is_root]

    if not branch_lengths:
        return {"max_edpl": 0.1, "max_pendant": 0.5}

    median_branch = statistics.median(branch_lengths)

    sorted_blens = sorted(branch_lengths)

    # 95th percentile of branch lengths. Values above this are likely in the realm
    # of LBA. Very broad/accepting approach to asses branches
    blen_95prc = sorted_blens[int(0.95 * (len(sorted_blens) - 1))]

    return median_branch * 2.0, max(blen_95prc * 2.0, median_branch * 5.0)


def classify_placement_quality(
        unc_smry: dict,
        max_pendant: float,
        max_edpl: float,
        min_agg_lwr: float = 0.80,
        max_entropy: float = 0.80) -> list:
    """
    Job is to simply pick out sequences with poor placement quality metrics.

    This is a pre-screening approach for the phylogenomic sistsers assessment.
    """

    place_qual, qual_info = None, None

    if unc_smry is None:
        place_qual = 'UNCERTAIN'
        qual_info = 'No-Placements'

    # Check if multiple placements span the proposed root
    elif unc_smry['LCA_Node'].is_root:
        place_qual = "UNCERTAIN"
        qual_info = "Placements-Span-Root"

    # Check LBA or possible novel-lineage
    elif unc_smry['Mean_Pendant'] > max_pendant:
        place_qual = 'UNCERTAIN'
        qual_info = 'Excessive-Branch-Length'

    # Checks phylogenetic spread (uncertain placements across the tree)
    elif unc_smry['EDPL'] > max_edpl:
        place_qual = 'UNCERTAIN'
        qual_info = 'High-EDPL-Spread'

    # Suffering from indecision, where placement is just terribly uncertain
    # Can be locally clustered, but still suffering from poor positioning
    elif unc_smry['Normalized_Entropy'] > max_entropy and unc_smry['Total_Placements'] > 2:
        place_qual = 'UNCERTAIN'
        qual_info = 'High-Placement-Entropy'

    # Cumulative Likelihood Weight Ratio check. Want to ensure the last common ancestor
    # is fairly certin
    elif unc_smry['LCA_Agg_LWR'] < min_agg_lwr:
        place_qual = 'UNCERTAIN'
        qual_info = 'Low-Aggregative-LWR'

    # Distinguish clean single-branch hits from local clade-level placements
    elif unc_smry['Total_Placements'] == 1 or unc_smry['Normalized_Entropy'] <= 0.30:
        place_qual = 'GOOD'
        qual_info = 'High-Confidence-Single-Branch'

    else:
        place_qual = 'GOOD'
        qual_info = 'High-Confidence-Clade-Consensus'

    unc_smry['Placement_Quality'] = place_qual
    unc_smry['Placement_Quality_Info'] = qual_info

    return unc_smry


def flatten_sister_taxonomy(query_txmy, unique_sister_txmy):
    # all_txnmy = set(query_txmy) | unique_sister_txmy
    query_txmy_levels = [[level] for level in query_txmy.split(';')]

    if len(query_txmy_levels) == 6:
        txmy_levels = ['Domain', 'Deep', 'Major', 'Minor', 'Shallow', 'Genus']

    else:
        return None

    prep_sister_txmy = [[level for level in txmy.split(";")] for txmy in unique_sister_txmy]

    sister_txmy_breakdown = {key: list(dict.fromkeys(column)) for key, column in zip(txmy_levels, zip(*prep_sister_txmy))}

    for k, v in sister_txmy_breakdown.items():
        if len(v) > 19:
            sister_txmy_breakdown[k] = 'too-many'
        else:
            sister_txmy_breakdown[k] = ';'.join(v)

    return sister_txmy_breakdown


def assess_sisters(query_txmy, all_sister_taxa, all_sister_seqs, unique_sister_txmy):
    """Checks start_node, then iterates upward using ete4's .ancestors() if needed."""

    if len(set(all_sister_seqs)) > 19:
        sisters_seq_str = 'too-many'

    else:
        sisters_seq_str = ';'.join(list(set(all_sister_seqs)))

    if len(set(all_sister_taxa)) > 19:
        sisters_taxa_str = 'too-many'

    else:
        sisters_taxa_str = ';'.join(list(set(all_sister_taxa)))

    sister_summary = flatten_sister_taxonomy(query_txmy, unique_sister_txmy)

    if not sister_summary:
        print('Error with taxonomic levels. Expecting six levels...')


    sister_summary['Total Sister Taxa'] = len(set(all_sister_taxa))
    sister_summary['Sister Taxa'] = sisters_taxa_str

    sister_summary['Total Sister Sequences'] = len(set(all_sister_seqs))
    sister_summary['Sister Sequences'] = sisters_seq_str

    return sister_summary


def prep_leaf_taxonomy(tree, query_taxon, delim: str = '_XX_') -> dict:
    leaf_taxonomy = nt.get_ncbi_txmy([query_taxon])

    for leaf in [i.name for i in tree.leaves()]:
        leaf_tmp_dict = nt.get_ncbi_txmy([leaf.partition(delim)[0]]) #nt.get_ncbi_txmy([leaf.partition(delim)[0]])

        if not leaf_tmp_dict:
            continue

        leaf_taxonomy.update(leaf_tmp_dict)
        leaf_taxonomy[leaf] = list(leaf_tmp_dict.values())[0]

    return leaf_taxonomy


def get_leaf_taxa(node, leaf_taxonomy, delim) -> List:
    taxa_list = []
    txnmy_list = []
    seq_list = []
    for i in node.leaves():
        taxa_list.append(i.name.partition(delim)[0])
        seq_list.append(i.name)
        txnmy_list.append(leaf_taxonomy[i.name])

    return [taxa_list, seq_list, set(txnmy_list)]


def assess_query_seqs(
        query_taxon: str,
        jplace_file: str,
        outdir: str,
        lwr_cutoff: float = 0.80,
        max_entropy: float = 0.8,
        delim: str = '_XX_'):

    tree, by_query = parse_jplace_file(jplace_file)

    leaf_taxonomy = prep_leaf_taxonomy(
                        tree,
                        query_taxon,
                        delim)

    max_edpl, max_pendant = tree_thresholds(tree)

    edpl_dict = load_gappa_edpl(
                    outdir,
                    jplace_file)

    diag_gf_summary = {}
    for k, v in by_query.items():
        query_taxon = k.partition(delim)[0]
        query_txmy = leaf_taxonomy[query_taxon]
        diag_gf_summary[k] = {'Diag_GF': k.rpartition(delim)[-1]}
        diag_gf_summary[k]['Max_EDPL'] = max_edpl
        diag_gf_summary[k]['Max_Pendant'] = max_pendant

        unc_smry = query_uncertainty_summary(
                        k,
                        tree,
                        v,
                        edpl_dict,
                        lwr_cutoff)

        query_place_smry = classify_placement_quality(
                                unc_smry,
                                max_pendant,
                                max_edpl,
                                lwr_cutoff,
                                max_entropy)

        diag_gf_summary[k].update(query_place_smry)

        for_assess_sis = [query_txmy] + get_leaf_taxa(query_place_smry['LCA_Node'], leaf_taxonomy, delim)

        query_txmy, all_sister_taxa, all_sister_seqs, unique_sister_txmy = for_assess_sis

        diag_gf_summary[k].update(assess_sisters(*for_assess_sis))

    return diag_gf_summary


def write_summary_table(
        tbl_outname: str,
        query_txmy_level: str,
        clade_level: str,
        eval_dict: dict
    ) -> None:

    with open(f'{tbl_outname}', 'w+') as w:
        w.write(f'Sister_Clade_Name,Sister_Counts,Clade_Type,Clade_Assesment\n')
        for k, v in eval_dict.items():
            if query_txmy_level == k:
                w.write(f'{k},{v},{clade_level},SAME\n')

            else:
                w.write(f'{k},{v},{clade_level},DIFF\n')



def merge_summary_tables(
        outdir: str,
        taxon_name: str,
        query_txmy: dict,
        tbl_prefix: str):

    all_tables = [pd.read_csv(tbl) for tbl in glob.glob(f'{tbl_prefix}*_Clade_Counts.csv')]

    merged_df = pd.concat(all_tables, axis = 0, ignore_index = True)

    merged_df['Query_Taxon'] = taxon_name
    merged_df['Query_Taxonomy'] = ';'.join(query_txmy[:4])

    merged_df = merged_df[['Query_Taxon','Query_Taxonomy','Sister_Clade_Name','Sister_Counts','Clade_Type','Clade_Assesment']]

    merged_df = merged_df.set_index('Query_Taxon')

    sorted_df = merged_df.sort_values(by = ['Clade_Type', 'Sister_Counts'], ascending = [True, False])

    lg_smry_tbl = f'{outdir}{taxon_name}.Phylogenomic.CompleteSummary.Clade_Counts.csv'

    sorted_df.to_csv(lg_smry_tbl)

    return lg_smry_tbl



def clean_summary_dict(summary_dict):
    return {k: v for k, v in summary_dict.items() if k not in ['too-many', 'non-monophyletic']}


def prep_summary_table(outdir, taxon_name, complete_sister_summary, eval_levels: list = ['Domain', 'Deep', 'Major', 'Minor']):
    query_txmy = list(nt.get_ncbi_txmy([taxon_name]).values())[0].split(';')
    smry_dict = {}
    fin_dict = {}

    for k, v in complete_sister_summary.items():
        if not v['Placement_Quality'] == 'GOOD':
            continue

        tmp_dict = {}

        for lvl in eval_levels:
            if ';' in v[lvl]:
                tmp_dict[lvl] = 'non-monophyletic'
            else:
                tmp_dict[lvl] = v[lvl]

        smry_dict[k] = tmp_dict

    eval_df = pd.DataFrame(smry_dict).T

    domain_smry = clean_summary_dict(eval_df.Domain.value_counts().to_dict())

    euk_eval_df = eval_df[eval_df['Domain'] == 'Eukaryota']

    deep_smry = clean_summary_dict(euk_eval_df.Deep.value_counts()[:10].to_dict())
    mjr_smry = clean_summary_dict(euk_eval_df.Major.value_counts()[:10].to_dict())
    mnr_smry = clean_summary_dict(euk_eval_df.Minor.value_counts()[:10].to_dict())

    write_summary_table(
        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic.Domain_Clade_Counts.csv',
        query_txmy[0],
        'Domain',
        domain_smry)

    write_summary_table(
        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic.Euk_Deep_Clade_Counts.csv',
        query_txmy[1],
        'Deep_Euk_Clade',
        deep_smry)

    write_summary_table(
        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic.Euk_Major_Clade_Counts.csv',
        query_txmy[2],
        'Major_Euk_Clade',
        mjr_smry)

    write_summary_table(
        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic.Euk_Minor_Clade_Counts.csv',
        query_txmy[3],
        'Minor_Euk_Clade',
        mnr_smry)


    fin_summary_tbl = merge_summary_tables(
                        outdir,
                        taxon_name,
                        query_txmy,
                        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic.')


def summarize_eval_phylo(
        start_time,
        backup_dir: str,
        outdir: str,
        jplace_dict,
        taxon_name: str,
        delim: str = '_XX_',
        lwr_cutoff: float = 0.8,
        max_entropy: float = 0.8):

    complete_sister_summary = {}
    out_csv  = f'{outdir}Phylogeny_Based/{taxon_name}.Complete_Sister_Summary.csv'

    fin_eval = 0
    total_diag = len(jplace_dict)

    for k, v in jplace_dict.items():
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Evaluated {fin_eval} of {total_diag} Diagnostic Phylogenies', end = '\r')
        diag_eval = assess_query_seqs(
                        taxon_name,
                        v,
                        outdir,
                        lwr_cutoff,
                        max_entropy,
                        delim)

        complete_sister_summary.update(diag_eval)

        fin_eval += 1

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Evaluated {fin_eval} of {total_diag} Diagnostic Phylogenies')

    sis_smry_df = pd.DataFrame(complete_sister_summary).T.drop(columns = ['LCA_Node'])

    sis_smry_df.index.name = 'Query_Sequence'

    sis_smry_df.to_csv(out_csv)

    prep_summary_table(outdir, taxon_name, complete_sister_summary)

        # complete_sister_summary.update()
