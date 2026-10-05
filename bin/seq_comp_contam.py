#!/usr/bin/env python3

"""
Analyses and plots include composition-based assessments of sequences to aid in
distinguishing between target and non-target (i.e., by-catch) sequence data.

This is hastily put together and there is a fair bit of redundancy... this will
be addressed eventually...

Part of the "by-catch" pipeline.

Dependencies:
- Python packages:
  + BioPython
  + Codon-Bias
  + Pandas
  + Scikit-Learn
"""

import itertools, sys

import codonbias as cb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from collections import defaultdict
from Bio import SeqIO

from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler


def calc_gc3_degen(seq: str):
    degen_codons = ['CT', 'GT','TC','CC','AC','GC','CG','GG']
    gc3_degen_str = ''.join([i[-1] for i in seq_kmer_counts(seq, 3, False) if i[:2] in degen_codons])
    return calc_gc_content(gc3_degen_str)


def calc_enc_gc3(fasta_file: str, genetic_code: int = 1, weighted: bool = False):
    """
    Weighted refers to the Sun et al. 2013 (MBE) paper! By default, their "robust"
    implementation of Wright's 1990 ENc calculation is on (not toggled)
    """
    enc_gc3_dict = {}
    if not weighted:
        enc = cb.scores.EffectiveNumberOfCodons(genetic_code = genetic_code, mean = 'unweighted')
    else:
        enc = cb.scores.EffectiveNumberOfCodons(genetic_code = genetic_code)
    for i in SeqIO.parse(fasta_file, 'fasta'):
        enc_gc3_dict[i.id] = {'GC3-Degen': calc_gc3_degen(f'{i.seq}')}
        enc_gc3_dict[i.id]['ENc'] = enc.get_score(f'{i.seq}')
    return enc_gc3_dict


def norm_kmer_counts(kmer_counts: list, kmers_to_eval: list):
    total_kmers = len(kmer_counts)
    norm_kmers = {tet: kmer_counts.count(tet)/total_kmers for tet in kmers_to_eval}
    return norm_kmers


def seq_kmer_counts(seq: str, kmer: int = 4, overlap = True) -> list:
    if overlap:
        return [seq[n:n + kmer] for n in range(len(seq)-kmer)]
    else:
        return [seq[n:n + kmer] for n in range(0, len(seq)-kmer, kmer)]

def calc_gc_content(seq: str):
    return float(f'{100 * ((seq.lower().count("g") + seq.lower().count("c")) / len(seq)):.2f}')


def calc_gc12_gc3(seq_kmers: list):
    gc12 = ''.join([i[:3] for i in seq_kmers])
    gc3 = ''.join([i[-1] for i in seq_kmers])
    gc123_dict = {'GC12': calc_gc_content(gc12), 'GC3': calc_gc_content(gc3)}
    return gc123_dict


def calc_nuc_freq(fasta_file: str) -> dict:
    fasta_tnf = {}
    fasta_gc123 = defaultdict(dict)
    gc123_kmer_list = [''.join(i) for i in itertools.product(['A', 'C', 'G', 'T'], repeat = 3)]
    tnf_kmer_list = [''.join(i) for i in itertools.product(['A', 'C', 'G', 'T'], repeat = 4)]
    for i in SeqIO.parse(fasta_file, 'fasta'):
        seq_kmers = seq_kmer_counts(f'{i.seq}', 3, False)
        fasta_gc123[i.id] = calc_gc12_gc3(seq_kmers)
        tnf_dict = norm_kmer_counts(seq_kmer_counts(f'{i.seq}', 4), tnf_kmer_list)
        fasta_tnf[i.id] = tnf_dict
    return fasta_tnf, fasta_gc123


def pca_tnf_scores(tnf_df):
    tnf_scaled = RobustScaler().fit_transform(tnf_df)
    pca = PCA().set_output(transform = "pandas")
    pca = pca.fit(tnf_scaled)
    tnf_pca = pca.transform(tnf_scaled)
    eigvals_dict = {f'PC{n+1}':{'Proportion of Variance':pca.explained_variance_ratio_[n],
                'Cumulative Proportion':np.cumsum(pca.explained_variance_ratio_)[n]}
                for n in range(len(pca.explained_variance_ratio_))}
    tnf_eigvals_df = pd.DataFrame(eigvals_dict)
    return tnf_pca, eigvals_dict


def store_composition_summaries(fasta_file: str, out_csv: str, genetic_code: int = 1, weighted: bool = False):
    tnf_dict, gc123_dict = calc_nuc_freq(fasta_file)
    enc_gc3_dict = calc_enc_gc3(fasta_file, genetic_code, weighted)

    tnf_df = pd.DataFrame(tnf_dict).T
    gc123_df = pd.DataFrame(gc123_dict).T
    enc_gc3_df = pd.DataFrame(enc_gc3_dict).T

    tnf_df.index.name = 'Sequence'
    gc123_df.index.name = 'Sequence'
    enc_gc3_df.index.name = 'Sequence'

    codon_df = pd.merge(gc123_df, enc_gc3_df, left_index = True, right_index = True)

    tnf_df.to_csv(f'{out_csv}TNF_Scores.csv')
    codon_df.to_csv(f'{out_csv}Codon_Bias.csv')

    return tnf_df, codon_df


def null_enc_gc3():
    null_gc3 = np.linspace(0.01, 0.99, 500)
    null_enc = 2 + null_gc3 + (29 / (null_gc3**2 + (1 - null_gc3)**2))

    return null_gc3, null_enc


def plot_composition_scores(fasta_file: str, tnf_pca, eigvals_dict, codon_df):
    out_png = f'{fasta_file.rpartition("/")[-1].rpartition(".")[0]}.SeqComposition.png'

    pca_title = f'PCA Plot for:\n{fasta_file.rpartition("/")[-1].partition(".")[0].replace("_"," ")}'
    gc123_title = pca_title.replace("PCA Plot ", "GC12 x GC3 Plot ")
    enc_title = pca_title.replace("PCA Plot ", "ENC Plot ")

    null_gc3, null_enc = null_enc_gc3()

    fig, ax = plt.subplots(1,3, figsize=(12, 5))

    sns.lineplot(x = null_gc3*100, y = null_enc, color = 'black', lw = 2, ax = ax[0])
    # sns.scatterplot(data = codon_df, x = 'GC3-Degen', y = 'ENc', alpha = 0.5, ax = ax[0])
    sns.kdeplot(data = codon_df, x = 'GC3-Degen', y = 'ENc', fill = True, cmap = 'flare', ax = ax[0])
    ax[0].set_xlabel(f'%GC3 4-fold Degenerate Codon Positions')
    ax[0].set_ylabel(f'ENc')
    ax[0].set_ylim(20, 65)
    ax[0].set_title(enc_title)

    sns.kdeplot(data = codon_df, x = 'GC3', y = 'GC12', fill = True, cmap = 'flare', ax = ax[1])
    ax[1].set_xlabel(f'%GC3')
    ax[1].set_ylabel(f'Average %GC12')
    ax[1].set_title(gc123_title)

    # sns.scatterplot(data = tnf_pca, x = 'pca0', y = 'pca1')
    sns.kdeplot(data = tnf_pca, x = 'pca0', y = 'pca1', fill = True, cmap = 'rocket_r', ax = ax[2])
    ax[2].set_xlabel(f'PC1 ({100*eigvals_dict['PC1']['Proportion of Variance']:.2f}%)')
    ax[2].set_ylabel(f'PC2 ({100*eigvals_dict['PC2']['Proportion of Variance']:.2f}%)')
    ax[2].set_title(pca_title)

    plt.tight_layout()
# plt.show()

    plt.savefig(out_png, dpi = 300)



def calc_comp_stats(fasta_file: str):
# Calculate the tetranucleotide frequencies
    out_csv = f'{fasta_file.rpartition("/")[-1].rpartition(".")[0]}.'

    tnf_df, codon_df = store_composition_summaries(fasta_file, out_csv)

    tnf_pca, eigvals_dict = pca_tnf_scores(tnf_df)

    plot_composition_scores(fasta_file, tnf_pca, eigvals_dict, codon_df)



if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
    except:
        print('Usage:\n\n python3 seq_compy_contam.py [FASTA-FILE]\n')
        sys.exit(1)

    calc_comp_stats(fasta_file)
