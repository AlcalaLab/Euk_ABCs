#!/usr/bin/env python3

"""
'Master' script to control the contamination pipeline!
"""

import glob, os, subprocess, sys, time

from bin import phylo_contam_args as pea
from bin import prep_query_contam as pqc
from bin import prep_phylo_tree_contam as ptc
from bin import eval_phylo_contam as epc

from datetime import timedelta
from pathlib import Path


def eval_args(args):
    if not pea.double_check_args(args):
        sys.exit(1)

    else:
        return args


def prep_query_data(
        start_time,
        backup_dir: str,
        outdir: str,
        fasta_file: str,
        taxon_name: str,
        prot_db: str,
        delim: str = '_XX_',
        approach: str = 'phylo',
        min_len: int = 300,
        threads: int = 4,
        gcode: int = 1
        ) -> None:

    print('#--------------- Preparing Query Data ---------------#')
    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Removing short transcripts')

    sfilt_fasta = pqc.remove_short_seqs(
                        fasta_file,
                        taxon_name,
                        outdir,
                        min_len)


    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Removing rRNA sequence by-catch')
    rrna_fasta = pqc.remove_rrna_seqs(
                        sfilt_fasta,
                        taxon_name,
                        backup_dir,
                        outdir,
                        threads)

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Assigning sequences to diagnostic gene families')
    orfs_by_gf = pqc.assign_diag_gene_families(
                        rrna_fasta,
                        taxon_name,
                        prot_db,
                        backup_dir,
                        outdir,
                        delim,
                        threads)

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Finalizing ORF calls')
    fin_prep_dir, query_seqs = pqc.finalize_query_orfs(
                                    prot_db,
                                    orfs_by_gf,
                                    outdir,
                                    delim,
                                    approach,
                                    gcode)

    try:
        Path(sfilt_fasta).unlink()
    except FileNotFoundError:
        pass

    try:
        Path(f'{sfilt_fasta}.fai').unlink()
    except FileNotFoundError:
        pass

    return fin_prep_dir, query_seqs


def update_phylo_data(
        start_time,
        backup_dir: str,
        outdir: str,
        taxon_name: str,
        unaln_dir: str,
        query_seqs: list,
        ref_aln_dir: str,
        ref_tree_dir: str,
        delim: str = '_XX_',
        threads: int = 4,
        blen_mode: str = 'median',
        eval_mode: str = 'short'
        ):

    print('\n#---- Preparing Query Alignments and Phylogenies ----#')
    updated_phylos = ptc.all_query_aln(
                        start_time,
                        outdir,
                        unaln_dir,
                        query_seqs,
                        ref_aln_dir,
                        ref_tree_dir,
                        taxon_name,
                        threads
                        )

    epc.summarize_eval_phylo(
            start_time,
            backup_dir,
            outdir,
            updated_phylos,
            taxon_name,
            delim,
            blen_mode,
            eval_mode
            )


def phylo_eval_contam(
        init_outdir: str,
        fasta_file: str,
        taxon_name: str,
        taxon_code: str,
        prot_db: str,
        ref_aln_dir: str,
        ref_tree_dir: str,
        delim: str = '_XX_',
        approach: str = 'phylo',
        min_len: int = 300,
        threads: int = 4,
        gcode: int = 1
    ) -> None:

    start_time = time.time()

    if taxon_code:
        taxon_name = f'{taxon_name}_{taxon_code}'

    outdir = f'{taxon_name}_Contam_Eval/'
    backup_dir = 'Phylo_Contam_Backup/'

    if init_outdir:
        backup_dir = f'{init_outdir}_Contam_Eval/Phylo_Contam_Backup/'
        outdir = f'{init_outdir.rstrip("/")}_Contam_Eval/{taxon_name}_Contam_Eval/'

    if glob.glob(f'{outdir}*.tar.gz'):
        print('\n#----------------------------------------------------#')
        print(f"[{timedelta(seconds = round(time.time()-start_time))}] {taxon_name}'s phylogeny-based assessment is already finished")
        print('#----------------------------------------------------#')

    else:
        print('\n#----------------------------------------------------#')
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Starting phylogeny-based assessment of {taxon_name}')
        print('#----------------------------------------------------#')

        unaln_dir, query_seqs = prep_query_data(
                                    start_time,
                                    backup_dir,
                                    outdir,
                                    fasta_file,
                                    taxon_name,
                                    prot_db,
                                    delim,
                                    approach,
                                    min_len,
                                    threads,
                                    gcode)

        update_phylo_data(
            start_time,
            backup_dir,
            outdir,
            taxon_name,
            unaln_dir,
            query_seqs,
            ref_aln_dir,
            ref_tree_dir,
            delim,
            threads
            )

        print('\n#----------------------------------------------------#')
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Finished phylogeny-based assessment of {taxon_name}')
        print('#----------------------------------------------------#')

        subprocess.run(
            f'tar -zcvf {outdir.rstrip("/")}.tar.gz -C {outdir} .',
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
            )

        os.system(f'rm -rf {outdir}')

        print('\n#----------------------------------------------------#')
        print(f"[{timedelta(seconds = round(time.time()-start_time))}]  Finished compressing {taxon_name}'s phylogeny-based assessment data")
        print('#----------------------------------------------------#')



if __name__ == '__main__':

    args = eval_args(pea.capture_args())

    print('NEED A FUNCTION TO CHECK THE FOLDER PATHS!!!! make sure that the ref-msa and ref-tree folders are there!')

    if not Path(args.db).is_file():
        print(f'Cannot find {args.db} in current folder')
        sys.exit()

    if args.phylo:
        if not Path(args.ref_msa).is_dir():
            print(f'Cannot find {args.ref_msa}')
            sys.exit()

        if not Path(args.ref_trees).is_dir():
            print(f'Cannot find {args.ref_trees}')
            sys.exit()

    if not args.out_dir:
        args.out_dir = ''

    if Path(args.input).is_file():
        if not args.taxon_name:
            args.taxon_name = args.input.rpartition("/")[-1].partition(".")[0]
        else:
            args.taxon_name = '_'.join(args.taxon_name)

        phylo_eval_contam(
            args.out_dir,
            args.input,
            args.taxon_name,
            args.taxon_code,
            args.db,
            args.ref_msa,
            args.ref_trees,
            '_XX_',
            args.approach,
            300,
            args.threads,
            1
            )

    elif Path(args.input).is_dir():

        query_fasta_files = glob.glob(f'{args.input}*fasta')
        print(len(query_fasta_files))

        if not query_fasta_files:
            query_fasta_files = glob.glob(f'{args.input}/*fas')

        elif not query_fasta_files:
            query_fasta_files = glob.glob(f'{args.input}/*fa')

        elif not query_fasta_files:
            query_fasta_files = glob.glob(f'{args.input}/*fna')

        elif not query_fasta_files:
            print(f'Cannot find FASTA files in the {args.input} folder...')
            sys.exit(1)

        else:
            pass

        for fasta_file in query_fasta_files:
            taxon_name = fasta_file.rpartition("/")[-1].partition(".")[0]
            phylo_eval_contam(
                args.out_dir,
                fasta_file,
                taxon_name,
                args.taxon_code,
                args.db,
                args.ref_msa,
                args.ref_trees,
                '_XX_',
                args.approach,
                300,
                args.threads,
                1
                )
    else:
        print(f'\nError: Cannot find {args.input}. Double check the given PATH\n')
        sys.exit(1)
