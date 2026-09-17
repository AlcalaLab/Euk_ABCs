#!/usr/bin/env python3
import glob, shutil, time
import numpy as np
import pandas as pd

from collections import defaultdict
from datetime import timedelta
from ete4 import Tree
from pathlib import Path

from bin import ncbi_taxonomy as nt



def eval_clade_sizes(
        tree,
        leaf_taxonomy: dict,
        taxon_delim: str = '_XX_') -> dict:

    major_clades = ['Bacteria', 'Archaea', 'Opisthokonta','Archaeplastida','Amoebozoa','Discoba','Metamonada','SAR','Orphan']

    clade_sizes = {i:[] for i in ['Prok']+major_clades[2:]}

    for node in tree.descendants("postorder"):
        mjr_c = []

        node_seqs = [i.name for i in node.leaves() if 'Unid' not in i.name]

        node_taxa = list(set([i.partition(taxon_delim)[0] for i in node_seqs]))

        for i in node.leaves():
            if 'Unid' in i.name:
                continue

            node_taxa.append(i.name.partition(taxon_delim)[0])

            node_txnmy = leaf_taxonomy[i.name]

            if node_txnmy[0] in ['Bacteria','Archaea']:
                mjr_c.append(node_txnmy[0])

            else:
                mjr_c.append(node_txnmy[1])

        node_taxa = list(set(node_taxa))

        if len(mjr_c) > 1:
            if mjr_c.count('Bacteria') + mjr_c.count('Archaea') == len(mjr_c):
                clade_sizes['Prok'].append((node, len(mjr_c), len(node_taxa), ';'.join(node_seqs)))

            else:
                for clade in major_clades[2:]:
                    if mjr_c.count(clade) == len(mjr_c):
                        taxa = list(set([]))

                        clade_sizes[clade].append((node, len(mjr_c), len(node_taxa), ';'.join(node_seqs)))

    return clade_sizes


def adjust_tree_root(
        tree,
        leaf_taxonomy: dict):
    """
    need to support custom clades!!!!
    """
    tree.set_midpoint_outgroup()

    clade_sizes = eval_clade_sizes(tree, leaf_taxonomy)

    for k, v in clade_sizes.items():
        if v:
            largest_clade = max(v, key = lambda x: x[2])
            if largest_clade[0]:
                tree.set_outgroup(largest_clade[0])
                break

        else:
            continue

    return tree


def simple_tree_walk(
        query_taxon,
        node,
        delim) -> list:
    tmp_leaves = [i.name.partition(delim)[0] for i in node.leaves()]
    if len(set(tmp_leaves)) == 1:
        return simple_tree_walk(query_taxon, node.parent, delim)
    else:
        return [i for i in tmp_leaves if query_taxon not in i]


def flatten_lineage_info(
        lineage_info: list) -> str:
    if len(set(lineage_info)) == 1:
        return lineage_info[0]
    else:
        return 'non-monophyletic'


def eval_sisters(
        sister_seqs: list,
        delim: str = '_XX_') -> str:
    node_sister_summary = {}
    sister_taxa = list(set([i.partition(delim)[0] for i in sister_seqs]))
    sister_taxonomy = nt.ncbi_taxonomy(sister_taxa)
    node_sister_summary['All_Deep'] = ';'.join(list(set([v[0] for v in sister_taxonomy.values()])))
    node_sister_summary['All_Major'] = ';'.join(list(set([v[1] for v in sister_taxonomy.values()])))
    node_sister_summary['All_Minor'] = ';'.join(list(set([v[2] for v in sister_taxonomy.values()])))
    node_sister_summary['All_Shallow'] = ';'.join(list(set([v[3] for v in sister_taxonomy.values()])))
    node_sister_summary['Deep'] = flatten_lineage_info([v[0] for v in sister_taxonomy.values()])
    node_sister_summary['Major'] = flatten_lineage_info([v[1] for v in sister_taxonomy.values()])
    node_sister_summary['Minor'] = flatten_lineage_info([v[2] for v in sister_taxonomy.values()])
    node_sister_summary['Shallow'] = flatten_lineage_info([v[3] for v in sister_taxonomy.values()])
    return node_sister_summary


def eval_blens(
        blen: float,
        threshold: float) -> str:

    if blen <= threshold:
        return 'Short'
    else:
        return 'Long'


def eval_phylo(
        tree_file: str,
        query_taxon: str,
        delim: str = '_XX_') -> dict:
    leaf_taxonomy = {}

    taxon_sister_summary = defaultdict(dict)

    # try just the branch lengths of JUST the tips vs the entire tree's branch lengths...
    tmp_tree = Tree(tree_file)

    # resolve any polytomies in the single gene tree
    tmp_tree.resolve_polytomy(descendants = True)

    all_leaves = [i for i in tmp_tree.leaves()]

    for leaf in [i.name for i in all_leaves]:
        leaf_tmp_dict = nt.ncbi_taxonomy([leaf.partition(delim)[0]])
        if not leaf_tmp_dict:
            continue
        leaf_taxonomy.update(leaf_tmp_dict)
        leaf_taxonomy[leaf] = list(leaf_tmp_dict.values())[0]

    # root the tree (Proks, then Opis, then "Excavata", etc...)
    tree = adjust_tree_root(tmp_tree, leaf_taxonomy)

    query_leaves = [i for i in all_leaves if query_taxon in i.name and 'XX_Ref_' not in i.name]

    tip_branch_lengths = [i.dist for i in all_leaves]

    median_short_blen = float(np.median(tip_branch_lengths))

    mean_short_blen = float(np.mean(tip_branch_lengths))

    q1, q3 = np.percentile(np.array(tip_branch_lengths), [25, 75])
    quartile_short_blen = float(q1)

    for node in query_leaves:
        # print(node.name)
        query_seq = node

        sister_seqs = simple_tree_walk(query_taxon, node, delim)

        node_sis_eval = eval_sisters(sister_seqs, delim)

        taxon_sister_summary[query_seq.name]['Tip_Branch_Length'] = query_seq.dist
        taxon_sister_summary[query_seq.name]['Median_Branch_Length_Threshold'] = median_short_blen
        taxon_sister_summary[query_seq.name]['Mean_Branch_Length_Threshold'] = mean_short_blen
        taxon_sister_summary[query_seq.name]['Quartile_Branch_Length_Threshold'] = quartile_short_blen

        taxon_sister_summary[query_seq.name]['Median_Branch_Length_Type'] \
            = eval_blens(query_seq.dist, median_short_blen)

        taxon_sister_summary[query_seq.name]['Mean_Branch_Length_Type'] \
            = eval_blens(query_seq.dist, mean_short_blen)

        taxon_sister_summary[query_seq.name]['Quartile_Branch_Length_Type'] \
            = eval_blens(query_seq.dist, quartile_short_blen)

        taxon_sister_summary[query_seq.name].update(node_sis_eval)

        if len(sister_seqs) < 20:
            taxon_sister_summary[query_seq.name]['sister-seqs'] = ';'.join(sister_seqs)

        else:
            taxon_sister_summary[query_seq.name]['sister-seqs'] = 'too-many-seqs'

    return taxon_sister_summary


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
                w.write(f'{k},{v},{clade_level.rstrip("_Clade")},SAME\n')

            else:
                w.write(f'{k},{v},{clade_level.rstrip("_Clade")},DIFF\n')

def summarize_eval_phylo(
        start_time,
        backup_dir: str,
        outdir: str,
        tree_file_list: list,
        taxon_name: str,
        delim: str = '_XX_',
        blen_mode: str = 'median',
        eval_mode: str = 'short'):

    complete_sister_summary = {}
    fin_eval = 0

    query_txmy = nt.ncbi_taxonomy([taxon_name])

    csv_backup_dir = f'{backup_dir}Phylo_Contam_SpreadSheets/'
    Path(csv_backup_dir).mkdir(exist_ok = True, parents = True)
    Path(f'{csv_backup_dir}Complete_Contamination_Summary/').mkdir(exist_ok = True, parents = True)
    # Path(f'{csv_backup_dir}Domain_Clade_Summary/').mkdir(exist_ok = True, parents = True)
    # Path(f'{csv_backup_dir}Major_Clade_Summary/').mkdir(exist_ok = True, parents = True)
    # Path(f'{csv_backup_dir}Minor_Clade_Summary/').mkdir(exist_ok = True, parents = True)
    # Path(f'{csv_backup_dir}Shallow_Clade_Summary/').mkdir(exist_ok = True, parents = True)
    Path(f'{csv_backup_dir}Overall_Clade_Summary/').mkdir(exist_ok = True, parents = True)

    big_table_out = f'{tree_file_list[0].rpartition('Updated_Phylogenies')[0]}{taxon_name}.Complete_Phylogenomic_Assessment.csv'

    for tree_file in tree_file_list:
        if Path(tree_file).is_file() != True:
            continue
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Evaluated {fin_eval} of {len(tree_file_list)} Diagnostic Phylogenies', end = '\r')

        tmp_dict = eval_phylo(
            tree_file,
            taxon_name,
            delim)

        if tmp_dict == None:
            print(tree_file)
            break

        for k, v in tmp_dict.items():
            complete_sister_summary[k] = v

        fin_eval += 1

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Evaluated {fin_eval} of {len(tree_file_list)} Diagnostic Phylogenies')

    eval_df = pd.DataFrame(complete_sister_summary).T

    if blen_mode == 'median':
        blen_eval = 'Median_'
    elif blen_mode == 'mean':
        blen_eval = 'Mean_'
    else:
        blen_eval = 'Quartile_'

    deep_cts = eval_df[eval_df[f'{blen_eval}Branch_Length_Type'] == 'Short'].Deep.value_counts().to_dict()

    all_euk_only = eval_df[eval_df.Deep == 'Eukaryota']
    all_euk_mjr_cts = all_euk_only.Major.value_counts()[:10].to_dict()
    all_euk_mnr_cts = all_euk_only.Minor.value_counts()[:10].to_dict()
    all_euk_slw_cts = all_euk_only.Shallow.value_counts()[:10].to_dict()

    euk_only =  eval_df[(eval_df.Deep == 'Eukaryota') & (eval_df[f'{blen_eval}Branch_Length_Type'] == 'Short')]
    euk_mjr_cts = euk_only.Major.value_counts()[:10].to_dict()
    euk_mnr_cts = euk_only.Minor.value_counts()[:10].to_dict()
    euk_slw_cts = euk_only.Shallow.value_counts()[:10].to_dict()

    write_summary_table(
        f'{outdir}{taxon_name}.Phylogenomic.Domain_Clade_Counts.csv',
        query_txmy[taxon_name][0],
        'Domain',
        deep_cts)

    write_summary_table(
        f'{outdir}{taxon_name}.Phylogenomic.Euk_Major_Clade_Counts.csv',
        query_txmy[taxon_name][1],
        'Major_Euk_Clade',
        euk_mjr_cts)

    write_summary_table(
        f'{outdir}{taxon_name}.Phylogenomic.Euk_Minor_Clade_Counts.csv',
        query_txmy[taxon_name][2],
        'Minor_Euk_Clade',
        euk_mnr_cts)

    write_summary_table(
        f'{outdir}{taxon_name}.Phylogenomic.Euk_Shallow_Clade_Counts.csv',
        query_txmy[taxon_name][-1],
        'Shallow_Euk_Clade',
        euk_slw_cts)

    write_summary_table(
        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic_LBs.Euk_Major_Clade_Counts.csv',
        query_txmy[taxon_name][1],
        'Major_Euk_Clade',
        all_euk_mjr_cts)

    write_summary_table(
        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic_LBs.Euk_Minor_Clade_Counts.csv',
        query_txmy[taxon_name][2],
        'Minor_Euk_Clade',
        all_euk_mnr_cts)

    write_summary_table(
        f'{outdir}Phylogeny_Based/{taxon_name}.Phylogenomic_LBs.Euk_Shallow_Clade_Counts.csv',
        query_txmy[taxon_name][-1],
        'Shallow_Euk_Clade',
        all_euk_slw_cts)

    eval_df.index.name = 'Query_Seq'

    eval_df.to_csv(big_table_out, index = True)

    fin_summary_tbl = merge_summary_tables(
                        outdir,
                        taxon_name,
                        query_txmy,
                        f'{outdir}{taxon_name}.Phylogenomic.')

    shutil.copy2(
        big_table_out,
        f'{csv_backup_dir}Complete_Contamination_Summary/')

    # shutil.copy2(
    #     f'{outdir}{taxon_name}.Phylogenomic.Domain_Clade_Counts.csv',
    #     f'{csv_backup_dir}Domain_Clade_Summary/')
    #
    # shutil.copy2(
    #     f'{outdir}{taxon_name}.Phylogenomic.Euk_Major_Clade_Counts.csv',
    #     f'{csv_backup_dir}Major_Clade_Summary/')
    #
    # shutil.copy2(
    #     f'{outdir}{taxon_name}.Phylogenomic.Euk_Minor_Clade_Counts.csv',
    #     f'{csv_backup_dir}Minor_Clade_Summary/')
    #
    # shutil.copy2(
    #     f'{outdir}{taxon_name}.Phylogenomic.Euk_Shallow_Clade_Counts.csv',
    #     f'{csv_backup_dir}Shallow_Clade_Summary/')

    shutil.copy2(
        fin_summary_tbl,
        f'{csv_backup_dir}Overall_Clade_Summary/')


def merge_summary_tables(
        outdir: str,
        taxon_name: str,
        query_txmy: dict,
        tbl_prefix: str
        ):
    all_tables = [pd.read_csv(tbl) for tbl in glob.glob(f'{tbl_prefix}*_Clade_Counts.csv')]

    merged_df = pd.concat(all_tables, axis = 0, ignore_index = True)

    merged_df['Query_Taxon'] = taxon_name
    merged_df['Query_Taxonomy'] = ';'.join(query_txmy[taxon_name][:4])

    merged_df = merged_df[['Query_Taxon','Query_Taxonomy','Sister_Clade_Name','Sister_Counts','Clade_Type','Clade_Assesment']]

    merged_df = merged_df.set_index('Query_Taxon')

    merged_df.sort_values(by = ['Clade_Type', 'Sister_Counts'], ascending = [True, False])

    lg_smry_tbl = f'{outdir}{taxon_name}.Phylogenomic.CompleteSummary.Clade_Counts.csv'

    merged_df.to_csv(lg_smry_tbl)

    return lg_smry_tbl
