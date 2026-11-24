#!/usr/bin/env python3

from ete4 import NCBITaxa
from pathlib import Path


ncbi = NCBITaxa()

def make_unclassified(label: str) -> str:
    if label.startswith('Unclassified_'):
        return label

    else:
        return f'Unclassified_{label}'


def make_fallback_lineage(include_species: bool = False, eukaryotic: bool = True, num_ranks: int = 8) -> str:
    if include_species:
        ranks_count = num_ranks

    else:
        ranks_count = num_ranks - 1

    if eukaryotic:
        levels = ['Eukaryota', 'Eukaryota_incertae_sedis'] + ['Unclassified_Eukaryota_incertae_sedis'] * ranks_count

    else:
        levels = ['Unclassified'] * (ranks_count + 2)

    return ';'.join(levels)


def resolve_to_taxid(taxon_name: str, include_species: bool = False, eukaryotic: bool = True) -> int | None:
    if include_species:
        query_taxon_name = taxon_name.replace('_',' ')
    else:
        query_taxon_name = ' '.join(taxon_name.replace('_',' ').split(' ')[:-1])

    taxid_dict = ncbi.get_name_translator([query_taxon_name])

    if not taxid_dict and 'Candidatus' in taxon_name:
        query_taxon_name = query_taxon_name.replace('Candidatus ','')
        taxid_dict = ncbi.get_name_translator([query_taxon_name])

    if query_taxon_name in taxid_dict and taxid_dict[query_taxon_name]:
        taxids = taxid_dict[query_taxon_name]

        if len(taxids) == 1:
            return taxids[0]

        if eukaryotic:
            for tid in taxids:
                lineage_names = set(ncbi.get_taxid_translator(ncbi.get_lineage(tid)).values())

                if 'Eukaryota' in lineage_names:
                    return tid

        return taxids[0]

    return None


def prep_deep_taxonomy(
        all_lineage_names: list,
        rank_to_names: dict,
        include_species: bool = False) -> list[str]:

    if include_species:
        ranks_to_extract = ['phylum', 'class', 'order', 'family', 'genus', 'species']

    else:
        ranks_to_extract = ['phylum', 'class', 'order', 'family', 'genus']

    valid_supergroups = {
        'Opisthokonta': 'Opisthokonta',
        'Amoebozoa': 'Amoebozoa',
        'Apusomonadida': 'Apusomonadida',
        'Breviatea': 'Breviatea',
        'Amorphea': 'Amorphea',
        'Sar': 'SAR',
        'Archaeplastida': 'Archaeplastida',
        'Glaucocystophyceae': 'Archaeplastida',
        'Viridiplantae': 'Archaeplastida',
        'Rhodelphida': 'Archaeplastida',
        'Rhodophyta': 'Archaeplastida',
        'Excavata': 'Excavata',
        'Discoba': 'Discoba',
        'Metamonada': 'Metamonada',
        'Cryptophyceae': 'Cryptista',
        'Cryptista': 'Cryptista',
        'Haptista': 'Haptista',
        'CRuMs': 'CRuMs',
        'Provora': 'Provora',
        'Hemimastigophora': 'Hemimastigophora',
        'Ancyromonadida': 'Ancyromonadida',
        'Malawimonadida': 'Malawimonadida',
        'Microheliella': 'Cryptista',
        'Palpitomonas': 'Cryptista'}


    valid_subgroups = {
        # SAR subgroups
        'Stramenopiles': 'Stramenopiles',
        'Alveolata': 'Alveolata',
        'Rhizaria': 'Rhizaria',

        # Amorphea subgroups
        'Metazoa': 'Metazoa',
        'Fungi': 'Fungi',
        'Choanoflagellata': 'Choanoflagellata',
        'Rotosphaerida': 'Rotosphaerida',
        'Filasterea': 'Filasterea',
        'Ichthyosporea': 'Ichthyosporea',
        'Pluriformea': 'Pluriformea',
        'Aphelidiomycota': 'Aphelidiomycota',
        'Microsporidia': 'Microsporidia',
        'Discosea': 'Discosea',
        'Evosea': 'Evosea',
        'Tubulinea': 'Tubulinea',

        # Archaeplastida subgroups
        'Chloroplastida': 'Chloroplastida',
        'Viridiplantae': 'Chloroplastida',
        'Rhodophyta': 'Rhodophyta',
        'Glaucocystophyceae': 'Glaucophyta',
        'Glaucophyta': 'Glaucophyta',
        'Rhodelphida': 'Rhodelphida',

        # Discoba / Metamonada / Excavata subgroups
        'Euglenozoa': 'Euglenozoa',
        'Heterolobosea': 'Heterolobosea',
        'Jakobida': 'Jakobida',
        'Parabasalia': 'Parabasalia',
        'Fornicata': 'Fornicata',
        'Preaxostyla': 'Preaxostyla',

        # Haptista / Cryptista subgroups
        'Haptophyta': 'Haptophyta',
        'Centroplasthelida': 'Centroplasthelida',
        'Cryptophyceae': 'Cryptophyceae',
        'Katablepharidophyta': 'Katablepharidophyta',

        # "Orphan" subgroups
        'Telonemia': 'Telonemia',
        'Nebulidia': 'Nebulidia',
        'Nibbleridia': 'Nibbleridia',
        'Ancyromonadidae': 'Ancyromonadidae',
        'Planomonadidae': 'Planomonadidae'}


    # Infers the domain
    if 'domain' in rank_to_names.keys():
        domain = rank_to_names['domain']

    elif eukaryotic and 'domain' not in rank_to_names.keys():
        domain = rank_to_names['domain']

    else:
        domain = 'Unclassified'


    if domain == 'Eukaryota':
        subgroup = None

        # Grab the corresponding supergroup name, otherwise return incertae sedis.
        # adds 'Eukaryota_' at the beginning for easier processing with Gappa
        supergroup = next(
            (valid_supergroups[lineage_name] for lineage_name in reversed(all_lineage_names) if lineage_name in valid_supergroups),
            'Eukaryota_incertae_sedis')

        if supergroup != 'Eukaryota_incertae_sedis':
            subgroup = next(
                (valid_subgroups[lineage_name] for lineage_name in reversed(all_lineage_names) if lineage_name in valid_subgroups),
                f'{supergroup}_incertae_sedis')
        else:
            subgroup = 'Eukaryota_incertae_sedis'

    elif domain in ['Bacteria', 'Archaea']:
        # add the 'Kingdom' level in place of supergroup for prokaryotes
        if 'kingdom' in rank_to_names.keys():
            supergroup = rank_to_names['kingdom']
            subgroup = rank_to_names['phylum']

    else:
        supergroup = domain
        subgroup = f'Unclassified_{domain}'

    return [domain, supergroup, subgroup]


def get_reduced_taxonomy(
        taxon_name: str,
        include_species: bool = False,
        eukaryotic: bool = True) -> str:

    if include_species:
        ranks_to_extract = ['phylum', 'class', 'genus', 'species']

    else:
        ranks_to_extract = ['phylum', 'class', 'genus']

    taxid = resolve_to_taxid(taxon_name, eukaryotic = eukaryotic)

    if taxid is None:
        return make_fallback_lineage(include_species, eukaryotic, len(ranks_to_extract))

    try:
        lineage_ids = ncbi.get_lineage(taxid)

    except ValueError:
        return make_fallback_lineage(include_species, eukaryotic)

    taxid_rank_names = ncbi.get_taxid_translator(lineage_ids)

    taxid_rank_levels = ncbi.get_rank(lineage_ids)

    all_lineage_names = [taxid_rank_names[tid] for tid in lineage_ids]

    rank_to_names = {rank: taxid_rank_names[tid] for tid, rank in taxid_rank_levels.items() if rank != 'no rank'}

    reduced_taxonomy = prep_deep_taxonomy(all_lineage_names, rank_to_names)

    if reduced_taxonomy[0] in ['Bacteria', 'Archaea']:
        ranks_to_extract.pop(0)
        ranks_to_extract.insert(1, 'order')

    last_known_rank = reduced_taxonomy[-1]

    for rank in ranks_to_extract:

        rank_name = rank_to_names.get(rank)
        if rank_name:
            reduced_taxonomy.append(rank_name)
            last_known_rank = rank_name
        else:
            reduced_taxonomy.append(make_unclassified(last_known_rank))

    return ";".join(reduced_taxonomy)


def get_ncbi_txmy(
        taxon_list: list,
        include_species: bool = False,
        eukaryotic: bool = True):

    txnmy_db = {}
    fixed_names = {
        'Brandtodinium_nutriculum': 'Zooxanthella_nutriculum',
        'Spironema_cf_multiciliatum': 'Spironema_multiciliatum',
        'Sulfonecta_cf_uniserialis': 'Sulfonecta_uniserialis'}

    if not Path(f'{Path.home()}/.local/share/ete/taxdump.tar.gz').is_file():
        print('oh-no, missing NCBI taxonomy database!')
        print('preparing NCBI taxonomy database')
        ncbi.update_taxonomy_database()

    for taxon_name in taxon_list:
        if taxon_name in fixed_names:
            taxon_name = fixed_names[taxon_name]

        euk = True

        if taxon_name.count("_") == 0:
            taxon_name += '_sp'

        if 'Candidatus' in taxon_name:
            euk = False

        txnmy_db[taxon_name] = get_reduced_taxonomy(
                                taxon_name,
                                include_species,
                                euk)

    return txnmy_db
