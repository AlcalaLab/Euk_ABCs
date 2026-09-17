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
- CD-HIT
- Barrnap
- DIAMOND
"""

import os, shutil, subprocess, sys, time

from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def remove_short_seqs(
        fasta_file: str,
        taxon_name: str,
        outdir: str,
        min_len: int = 300
        ) -> str:

    sfilt_seqs = []
    renaming_dict = {}
    txp_num = 1

    outdir_shorts = f'{outdir}Data_Preparation/Size_Filtered/'
    out_sfilt_fasta = f'{outdir_shorts}{taxon_name}.{min_len}nt.fasta'
    out_sfilt_csv = f'{out_sfilt_fasta.replace(".fasta",".UpdatedNames.csv")}'

    Path(f'{outdir_shorts}').mkdir(exist_ok = True, parents = True)

    for i in SeqIO.parse(fasta_file, 'fasta'):
        if len(i.seq) < min_len:
            continue

        else:
            base_name = f'{taxon_name}_XX_Transcript_{txp_num}_Len_{len(i.seq)}'

            if 'cov_' in i.id:
                base_name += f'_Kcov_{float(i.id.partition("_cov_")[-1].partition("_")[0]):.2f}'

            renaming_dict[i.id] = base_name
            i.id = base_name
            i.description = ''
            i.name = ''

            sfilt_seqs.append(i)

            txp_num += 1

    SeqIO.write(sfilt_seqs, out_sfilt_fasta, 'fasta')

    with open(out_sfilt_csv, 'w+') as w:
        w.write('Initial_Name,Updated_Name\n')
        for k, v in renaming_dict.items():
            w.write(f'{k},{v}\n')

    return out_sfilt_fasta


def run_barrnap(
        fasta_file: str,
        taxon_name: str,
        outdir: str,
        threads: int = 4
    ) -> list:

    out_rrna_fastas = []

    outdir_rrna = f'{outdir}Data_Preparation/rRNA_Filtered/'

    out_bnap_base_fasta = f'{outdir_rrna}{taxon_name}.'

    kgdm = {'Eukaryotic':'euk','Bacterial':'bac','Mitochondrial':'mito','Archaeal':'arc'}

    for k, v in kgdm.items():
        rrna_fasta = f'{out_bnap_base_fasta}{k}_rRNA.fasta'
        bnap_cmd = f'barrnap --quiet ' \
            f'--threads {threads} ' \
            f'--kingdom {v} ' \
            f'--outseq {rrna_fasta} ' \
            f'{fasta_file}'

        out_rrna_fastas.append(rrna_fasta)

        subprocess.run(
            bnap_cmd,
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

    return out_rrna_fastas


def remove_rrna_seqs(
        fasta_file: str,
        taxon_name: str,
        backup_dir: str,
        outdir: str,
        threads: int = 4
        ) -> str:

    rrna_seqs = []

    outdir_rrna = f'{outdir}Data_Preparation/rRNA_Filtered/'

    out_rfilt_fasta = f'{outdir_rrna}{taxon_name}.rRNA_Filtered.fasta'

    tmp_rrna_fasta = f'{outdir_rrna}{taxon_name}.rRNA_Sequences.TMP.fasta'
    all_rrna_fasta = tmp_rrna_fasta.replace(".TMP.fasta",".fasta")

    Path(outdir_rrna).mkdir(exist_ok = True, parents = True)

    if Path(out_rfilt_fasta).is_file():
        if os.stat(out_rfilt_fasta).st_size != 0:
            return out_rfilt_fasta

    all_rrna_fastas = run_barrnap(
                        fasta_file,
                        taxon_name,
                        outdir,
                        threads)

    time.sleep(1)

    for f in all_rrna_fastas:
        rrna_seqs += [i for i in SeqIO.parse(f,'fasta')]

    SeqIO.write(rrna_seqs, tmp_rrna_fasta, 'fasta')

    time.sleep(1)

    rrna_seqs_to_toss = list(set([i.id.split("::")[1].rpartition(":")[0] for i in rrna_seqs]))

    cdhit_cmd = 'cd-hit-est -G 0 ' \
        '-c 0.99 ' \
        '-aS 1.0 ' \
        '-aL .0005 ' \
        f'-T {threads} ' \
        f'-i {tmp_rrna_fasta} ' \
        f'-o {all_rrna_fasta}'

    subprocess.run(
        cdhit_cmd,
        shell = True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
        )

    for f in all_rrna_fastas:
        subprocess.run(
            f'rm {tmp_rrna_fasta}',
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

    rrna_filt_seqs = [i for i in SeqIO.parse(fasta_file,'fasta') if i.id not in rrna_seqs_to_toss]

    SeqIO.write(rrna_filt_seqs, out_rfilt_fasta, 'fasta')

    rrna_backup_dir = f'{backup_dir}Phylo_Contam_FASTAs/rRNA_ByCatch_FASTAs/'

    Path(rrna_backup_dir).mkdir(exist_ok = True, parents = True)
    shutil.copy2(all_rrna_fasta, rrna_backup_dir)

    return out_rfilt_fasta


def extract_orfs(
        backup_dir: str,
        fasta_file: str,
        out_orf_tsv: str,
        delim: str = '_XX_',
        threads: int = 4
        ) -> str:

    og_assignments = []
    orfs_by_gf = defaultdict(list)

    out_orf_fasta = f'{out_orf_tsv.rpartition(".BLASTX")[0]}.Diag_ORFs.fasta'

    if Path(out_orf_fasta).is_file():
        if os.stat(out_orf_fasta).st_size != 0:
            for i in SeqIO.parse(out_orf_fasta,'fasta'):
                orfs_by_gf[i.id.rpartition(delim)[-1]].append(i)
            return orfs_by_gf

    init_seqs = {i.id:i.seq for i in SeqIO.parse(fasta_file, 'fasta')}

    for line in open(out_orf_tsv).readlines():
        qorf = ''

        qseqid, og, qstart, qend, qframe = line.split('\t')

        tmp_seq = init_seqs[qseqid]

        if int(qframe) > 0:
            qorf = tmp_seq[int(qstart)-1:int(qend)]
        else:
            qorf = tmp_seq[int(qend)-1:int(qstart)].reverse_complement()

        if qorf:
            seq_rec = SeqRecord(
                Seq(qorf),
                id = f'{qseqid}{delim}{og.rpartition(delim)[-1]}',
                description = '',
                name = ''
                )
            og_assignments.append(seq_rec)
            orfs_by_gf[og.rpartition(delim)[-1]].append(seq_rec)

    SeqIO.write(og_assignments, out_orf_fasta, 'fasta')

    orf_backup_dir = f'{backup_dir}Phylo_Contam_FASTAs/Diag_ORFs_FASTAs/'

    Path(orf_backup_dir).mkdir(exist_ok = True)
    shutil.copy2(out_orf_fasta, orf_backup_dir)

    return orfs_by_gf


def assign_diag_gene_families(
        fasta_file: str,
        taxon_name: str,
        prot_db: str,
        backup_dir: str,
        outdir: str,
        delim: str = '_XX_',
        threads: int = 4
        ) -> dict:

    outdir_orfs = f'{outdir}Data_Preparation/ORF_Calling/'

    out_orf_tsv = f'{outdir_orfs}{taxon_name}.Diag_Gene_Families.BLASTX.tsv'

    Path(outdir_orfs).mkdir(exist_ok = True, parents = True)

    if Path(out_orf_tsv).is_file():
        if os.stat(out_orf_tsv).st_size != 0:
            return extract_orfs(
                        backup_dir,
                        fasta_file,
                        out_orf_tsv,
                        delim,
                        threads
                        )

    dmnd_cmd = 'diamond blastx -k 1 -e 1e-10 --subject-cover 30 --very-sensitive ' \
        f'--threads {threads} ' \
        f'-q {fasta_file} ' \
        f'-d {prot_db} ' \
        f'-o {out_orf_tsv} ' \
        f'-f 6 qseqid sseqid qstart qend qframe'

    subprocess.run(
        dmnd_cmd,
        shell = True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
        )

    return extract_orfs(
            backup_dir,
            fasta_file,
            out_orf_tsv,
            delim,
            threads
            )


def prep_unaln_orfs(
        orfs_by_gf: dict,
        outdir: str,
        gcode: int = 1
    ) -> None:

    all_query_seqs = []

    for k, v in orfs_by_gf.items():
        tmp = []
        for i in v:
            tmp_seq = i.seq.translate(gcode)
            i.seq = tmp_seq
            tmp.append(i)
            all_query_seqs.append(i.id)

        SeqIO.write(tmp, f'{outdir}{k}.Query_AAs.fasta', 'fasta')

    return all_query_seqs


# def prep_af_orfs(
#         prot_db: str,
#         orfs_by_gf: dict,
#         outdir: str,
#         delim: str = '_XX_',
#         gcode: int = 1
#     ) -> None:
#
#     gf_dict = defaultdict(list)
#
#     for k, v in orfs_by_gf.items():
#         for i in v:
#             tmp_seq = i.seq.translate(gcode)
#             i.seq = tmp_seq
#             gf_dict[k].append(i)
#
#     for i in SeqIO.parse(prot_db,'fasta'):
#         gf = i.id.rpartition(delim)[-1]
#         # print(gf)
#         if gf in gf_dict.keys():
#             gf_dict[i.id.rpartition(delim)[-1]].append(i)
#
#     for k, v in gf_dict.items():
#         SeqIO.write(v, f'{outdir}{k}.AlignFree_GFs.fasta', 'fasta')


def finalize_query_orfs(
        prot_db: str,
        orfs_by_gf: dict,
        outdir: str,
        delim: str = '_XX_',
        approach: str = 'phylo',
        gcode: int = 1) -> str:

    if approach == 'phylo':
        prepped_outdir = f'{outdir}Phylogeny_Based/UnAligned_Query_GFs/'

    else:
        prepped_outdir = f'{outdir}Alignment_Free_Based/Query_GFs/'

    Path(prepped_outdir).mkdir(exist_ok = True, parents = True)

    if approach == 'phylo':
        query_seqs = prep_unaln_orfs(
            orfs_by_gf,
            prepped_outdir,
            gcode
            )
    else:
        prep_af_orfs(
            prot_db,
            orfs_by_gf,
            prepped_outdir,
            delim,
            gcode
            )

    return prepped_outdir, query_seqs
