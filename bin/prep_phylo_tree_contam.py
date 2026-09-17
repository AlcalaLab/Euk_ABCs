#!/usr/bin/env python3

"""
All the steps to prepare the data to assess for contamination.

Includes:
- removal of short transcripts
- removal and clustering of rRNA sequences
- assignment to diagnostic gene families
- extraction of ORFs from transcripts

Needs:
- BioPython
- MAFFT
- EPA-ng
- GAPPA
"""

import glob, os, subprocess, sys, time

from datetime import timedelta
from pathlib import Path

from Bio import SeqIO


def query_to_aln(
        fasta_file: str,
        ref_aln_fasta: str,
        outdir: str,
        query_seqs: list,
        threads: int = 4
        ) -> dict:

    query_aln_dir = f'{outdir}Phylogeny_Based/Updated_Alignments/'

    Path(query_aln_dir).mkdir(exist_ok = True, parents = True)

    out_aln_fasta = f'{query_aln_dir}{fasta_file.rpartition("/")[-1].rpartition(".")[0]}.Ref_Aln.fasta'

    mafft_cmd = f'mafft --quiet ' \
                f'--thread {threads} ' \
                '--auto ' \
                f'--add {fasta_file} ' \
                f'--keeplength {ref_aln_fasta} ' \
                f'> {out_aln_fasta}'

    if not (Path(out_aln_fasta).is_file() and Path(out_aln_fasta).stat().st_size > 0):
        subprocess.run(
            mafft_cmd,
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

    query_seqs = [i for i in SeqIO.parse(out_aln_fasta, 'fasta') if i.id in query_seqs]

    SeqIO.write(query_seqs, out_aln_fasta, 'fasta')

    return out_aln_fasta


def seqs_to_ref_tree(
        outdir: str,
        query_aln: str,
        ref_msa: str,
        ref_tree: str,
        threads: int = 4
    ) -> str:

    updated_tree_dir = f'{outdir}Phylogeny_Based/Updated_Phylogenies/'

    jplace_file = f'{updated_tree_dir}{query_aln.rpartition("/")[-1].rpartition(".")[0]}.EPAng.jplace'
    out_jplace = jplace_file.replace('.fasta_Ref.Aln.','.Ref_Aln.')
    out_newick = out_jplace.replace(".jplace",".newick")

    Path(updated_tree_dir).mkdir(exist_ok = True, parents = True)
    # apples_cmd = 'run_apples.py -p ' \
    #             f'-T {threads} ' \
    #             f'-t {ref_tree} ' \
    #             f'-s {ref_msa} ' \
    #             f'-x {query_aln} ' \
    #             f'-o {out_jplace}'

    epang_cmd = 'epa-ng -m LG ' \
                f'-T {threads} ' \
                f'-t {ref_tree} ' \
                f'-s {ref_msa} ' \
                f'-q {query_aln} ' \
                f'-w {updated_tree_dir}'

    gappa_cmd = 'gappa ' \
                'examine ' \
                'graft ' \
                '--fully-resolve ' \
                f'--threads {threads} ' \
                f'--jplace-path {out_jplace} ' \
                f'--out-dir {updated_tree_dir}'


    if not (Path(out_newick).is_file() and Path(out_newick).stat().st_size > 0):

        subprocess.run(
            epang_cmd,
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

        subprocess.run(
            f'mv {updated_tree_dir}epa_result.jplace {out_jplace}',
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

        os.system(f'rm {updated_tree_dir}epa_*log')

        subprocess.run(
            gappa_cmd,
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

        os.system(f'rm {updated_tree_dir}*jplace')

    return out_newick


def all_query_aln(
    start_time,
    outdir: str,
    outdir_unaln: str,
    query_seqs: list,
    ref_aln_dir: str,
    ref_tree_dir: str,
    query_taxon: str,
    threads: int = 4
    ) -> None:

    diag_db = {}
    updated_phylos = []
    with_updated_aln = 0

    for f in glob.glob(f'{outdir_unaln}*.fasta'):
        diag_gf = f.rpartition("/")[-1].partition(".")[0]
        diag_db[diag_gf] = [f]

    for f in glob.glob(f'{ref_aln_dir}*fasta'):
        # print(f)
        diag_gf = f.rpartition("/")[-1].partition(".")[0]
        if diag_gf in diag_db.keys():
            diag_db[diag_gf].append(f)

    for f in glob.glob(f'{ref_tree_dir}*nwk'):
        # print(f)
        diag_gf = f.rpartition("/")[-1].partition(".")[0]
        if diag_gf in diag_db.keys():
            diag_db[diag_gf].append(f)

    for k, v in diag_db.items():
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Updated {with_updated_aln} of {len(diag_db)} Reference Alignments', end = '\r')
        query_aln = query_to_aln(
                        v[0],
                        v[1],
                        outdir,
                        query_seqs,
                        threads
                        )
        # diag_db[k] += [i for i in SeqIO.parse(query_aln, 'fasta')]
        diag_db[k].append(query_aln)
        with_updated_aln += 1

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Updated {with_updated_aln} of {len(diag_db)} Reference Alignments')

    for k, v in diag_db.items():
        # print(k, v)
        if len(v) < 4:
            print(k, v)
            # print(v)
            sys.exit()
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Updated {len(updated_phylos)} of {len(diag_db)} Reference Phylogenies', end = '\r')
        updated_phylos.append(
                seqs_to_ref_tree(
                outdir,
                v[-1],
                v[1],
                v[2],
                threads)
                )

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Updated {len(updated_phylos)} of {len(diag_db)} Reference Phylogenies')

    return updated_phylos
