#!/usr/bin/env python3

from ete4 import NCBITaxa
from pathlib import Path

useful_ranks = ['domain','clade','phylum','family']

def ncbi_taxonomy(
        taxon_list: list,
        update_ncbi: bool = False) -> dict:

    txnmy_db = {}
    ncbi = NCBITaxa()

    if not Path(f'{Path.home()}/.local/share/ete/taxdump.tar.gz').is_file():
        print('oh-no, missing NCBI taxonomy database!')
        print('preparing NCBI taxonomy database')
        ncbi.update_taxonomy_database()

    for taxon in taxon_list:
        if 'Unid' in taxon:
            continue

        if 'Candidatus' in taxon:
            genus = taxon.split("_")[1]

        else:
            genus = taxon.partition("_")[0].replace("Prasinophyte","Volvox").replace("Brandtodinium","Zooxanthella").replace('Rhodophyte',"Galdieria")
            genus = genus.replace("Cercomonad","Cercozoa").replace('Cryptophyte','Cryptophyta').replace("Prymnesium","Volvox")
        # print(genus)
        clear_taxon = True

        try:
            taxid = list(ncbi.get_name_translator([genus]).values())[0][0]

        except IndexError:
            clear_taxon = False
            txnmy_db[taxon] = ['Eukaryota','Incertae sedis','Incertae sedis', 'Incertae sedis']

        if clear_taxon:
            tax_lineage = ncbi.get_lineage(taxid)
            full_lineage_names = ncbi.get_taxid_translator(tax_lineage)
            full_taxonomy = [full_lineage_names[taxid] for taxid in tax_lineage]
            txnmy_db[taxon] = reduce_taxonomy(full_taxonomy, genus)

    return txnmy_db


def reduce_taxonomy(
        full_taxonomy: list,
        taxon_name: str) -> list:

    if 'Opisthokonta' in full_taxonomy:
        if 'Metazoa' in full_taxonomy:
            reduced_taxonomy = full_taxonomy[2:5]+[full_taxonomy[7]]
        elif 'Fungi' in full_taxonomy:
            reduced_taxonomy = full_taxonomy[2:5]+[full_taxonomy[6]]
        else:
            reduced_taxonomy = full_taxonomy[2:6]

    elif 'Viridiplantae' in full_taxonomy:
        reduced_taxonomy = [full_taxonomy[2],'Archaeplastida', 'Chloroplastida', full_taxonomy[6]]

    elif 'Rhodophyta' in full_taxonomy:
        if len(full_taxonomy[2:]) < 4:
            reduced_taxonomy = ['Eukaryota', 'Archaeplastida', 'Rhodophyta incertae sedis', taxon_name]
        else:
            reduced_taxonomy = [full_taxonomy[2],'Archaeplastida', 'Rhodophyta', full_taxonomy[6]]

    elif 'Rhodelphea' in full_taxonomy:
        reduced_taxonomy = [full_taxonomy[2], 'Archaeplastida', full_taxonomy[3], full_taxonomy[4]]

    elif 'Cryptophyceae' in full_taxonomy:
        if len(full_taxonomy[2:]) < 4:
            reduced_taxonomy = ['Eukaryota','Cryptista', 'Cryptista incertae sedis', taxon_name]
        else:
            reduced_taxonomy = ['Eukaryota','Cryptista', full_taxonomy[4], full_taxonomy[5]]

    elif 'Sar' in full_taxonomy:
        if ('Ochrophyta' or 'Bigrya') in full_taxonomy[5]:
            reduced_taxonomy = ['Eukaryota', 'SAR', full_taxonomy[4], full_taxonomy[6]]
        else:
            reduced_taxonomy = ['Eukaryota', 'SAR', full_taxonomy[4], full_taxonomy[5]]

    elif len(full_taxonomy[2:]) < 4:
        reduced_taxonomy = full_taxonomy[2:-1]+[f'{full_taxonomy[2]} incertae sedis']+[full_taxonomy[-1]]

    else:
        reduced_taxonomy = full_taxonomy[2:6]

    return reduced_taxonomy
