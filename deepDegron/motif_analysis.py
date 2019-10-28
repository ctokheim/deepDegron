"""
File: motif_analysis.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Perform a motif analysis of deepDegron results
"""

import degron_pred
import numpy as np
import scipy.stats as stats
import pandas as pd
import argparse
import os
import pickle

import sklearn
import sklearn.metrics as metrics
from kneed import KneeLocator

import pvalue


def parse_arguments():
    info = 'Perform a motif analysis'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Results of deepDegron prediction')
    parser.add_argument('-t', '--type',
                        type=str, required=True,
                        help='Type of deepDegron prediction')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='motif enrichment ouput directory')
    args = parser.parse_args()
    return vars(args)


def compare_motif_cterm(peptide, motif):
    """C-term position specific motif match."""
    for i in range(len(motif)):
        if peptide[-(i+1)]!=motif[-(i+1)] and motif[-(i+1)]!='x':
            return 0
    return 1


def compare_motif_cterm_nopos(peptide, motif):
    """C-term not position specific motif match."""
    motif_ct = peptide[-6:].count(motif)
    return motif_ct


def compare_motif_nterm(peptide, motif):
    """N-term position specific motif match."""
    for i in range(len(motif)):
        if peptide[i]!=motif[i] and motif[i]!='x':
            return 0
    return 1


def compare_motif_nterm_nopos(peptide, motif):
    """C-term not position specific motif match."""
    motif_ct = peptide[:6].count(motif)
    return motif_ct


def motif_count(mystring, motif_list, motif_type='cterm'):
    """Count position specific motifs"""
    mysum = 0
    for motif in motif_list:
        if motif_type == 'cterm':
            mysum += compare_motif_cterm(mystring, motif)
        else:
            mysum += compare_motif_nterm(mystring, motif)
    return mysum


def motif_count_nopos(mystring, motif_list, motif_type='cterm'):
    mysum = 0
    for motif in motif_list:
        if motif_type == 'cterm':
            mysum += compare_motif_cterm_nopos(mystring, motif)
        else:
            mysum += compare_motif_nterm_nopos(mystring, motif)
    return mysum


def calc_top_degron_threshold(result):
    """See how many top degron potential sequences need to be dropped before
    the bag of amino acids and position specific models provide similar results."""
    step = 20
    possible_thresh = np.arange(step, 8000, step)

    # get the delta auc for every threshold
    output_list = []
    for thresh in possible_thresh:
        tmp = result.iloc[thresh:-thresh]
        pos_auc = metrics.roc_auc_score(tmp['y'], tmp['sequence position specific'])
        bag_auc = metrics.roc_auc_score(tmp['y'], tmp['bag of words'])
        delta_auc = pos_auc - bag_auc
        output_list.append([thresh, delta_auc])
        #if delta_auc < 0.01: return thresh-step

    # figure out the knee
    result_df = pd.DataFrame(output_list, columns=['threshold', 'delta auc'])
    knee_obj = KneeLocator(result_df['threshold'], result_df['delta auc'],
                           curve='convex', direction='decreasing')

    return knee_obj.knee


def nterm_dimotif_enrichment(df, thresh):
    aa_vocab = degron_pred.vocab
    output_list = []
    base_str = 'xxxxxx'
    for i in range(1,6):
        for j in range(i+1,7):
            for aa1 in aa_vocab:
                for aa2 in aa_vocab:
                    tmp = list(base_str)
                    tmp[-i] = aa1
                    tmp[-j] = aa2
                    pattern = ''.join(tmp)
                    testing = df['Peptide amino acid sequence'].apply(motif_count, args=([pattern], 'nterm'))
                    top_ct = testing.iloc[:thresh].sum()
                    #bottom_ct = testing.iloc[thresh:-thresh].sum() + 1
                    bottom_ct = testing.sum()
                    top_n = thresh
                    #bot_n = len(testing)-thresh*2+1
                    bot_n = len(testing)
                    pval = stats.binom_test(top_ct, n=top_n, p=bottom_ct/bot_n, alternative='greater')
                    #output_list.append([pattern, i, j, aa1, aa2, top_ct, bottom_ct, thresh, len(testing)-2*thresh, pval])
                    output_list.append([pattern, i, j, aa1, aa2, top_ct, bottom_ct, thresh, len(testing), pval])

    # compile results
    mycols = ['motif', 'pos1', 'pos2', 'aa1', 'aa2', 'top ct', 'bot ct', 'top total', 'bot total', 'pvalue']
    result = pd.DataFrame(output_list, columns=mycols)
    result['log(OR)'] = np.log2((result['top ct'] / (result['top total']-result['top ct'])) / (result['bot ct'] / (result['bot total']-result['bot ct'])))
    result['adjusted pvalue'] = result['pvalue'] * len(result)
    result['qvalue'] = pvalue.bh_fdr(result['pvalue'])
    result.sort_values('pvalue', inplace=True)

    return result


def nterm_trimotif_enrichment(df, thresh):
    aa_vocab = degron_pred.vocab
    # get background probs for each aa
    aa_bg = {}
    #bot_n = len(df)-thresh*2+1
    bot_n = len(df)
    for aa in aa_vocab:
        tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([aa], 'nterm'))
        #tmp_prob = (tmp.iloc[thresh:-thresh].sum() + 1) / (bot_n*6)
        tmp_prob = (tmp.sum()) / (bot_n*6)
        aa_bg[aa] = tmp_prob

    output_list = []
    for aa1 in aa_vocab:
        for aa2 in aa_vocab:
            for aa3 in aa_vocab:
                pattern = ''.join([aa1, aa2, aa3])
                # figure out full count
                testing = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern], 'nterm'))
                top_ct = testing.iloc[:thresh].sum()
                top_n = thresh * 4
                # get baseline prob
                bottom_prob_1 = aa_bg[pattern[0]]
                tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern[:2]], 'nterm'))
                #bottom_prob_2_di = (tmp.iloc[thresh:-thresh].sum() + 1) / (bot_n*5)
                bottom_prob_2_di = (tmp.sum()) / (bot_n*5)
                bottom_prob_2 = bottom_prob_2_di / bottom_prob_1
                bottom_prob_2_aa = aa_bg[pattern[1]]
                tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern[1:3]], 'nterm'))
                #bottom_prob_3_di = (tmp.iloc[thresh:-thresh].sum() + 1) / (bot_n*5)
                bottom_prob_3_di = (tmp.sum()) / (bot_n*5)
                bottom_prob_3 = bottom_prob_3_di / bottom_prob_2_aa
                bottom_prob = bottom_prob_1 * bottom_prob_2 * bottom_prob_3

                # measure significance
                pval = stats.binom_test(top_ct, n=top_n, p=bottom_prob, alternative='greater')
                #output_list.append([pattern, top_ct, bottom_prob, thresh, len(testing)-2*thresh, pval])
                output_list.append([pattern, top_ct, bottom_prob, thresh, bot_n, pval])

    # compile results
    mycols = ['motif', 'top ct', 'background p', 'top total', 'bot total', 'pvalue']
    result = pd.DataFrame(output_list, columns=mycols)
    result['log(OR)'] = np.log2((result['top ct'] / (result['top total']-result['top ct'])) / (result['background p'] / (1 - result['background p'])))
    result['adjusted pvalue'] = result['pvalue'] * len(result)
    result.loc[result['adjusted pvalue']>1, 'adjusted pvalue'] = 1
    result['qvalue'] = pvalue.bh_fdr(result['pvalue'])
    result.sort_values('pvalue', inplace=True)

    return result


def cterm_dimotif_enrichment(df, thresh):
    aa_vocab = degron_pred.vocab
    output_list = []
    base_str = 'xxxxxx'
    for i in range(1,6):
        for j in range(i+1,7):
            for aa1 in aa_vocab:
                for aa2 in aa_vocab:
                    tmp = list(base_str)
                    tmp[-i] = aa1
                    tmp[-j] = aa2
                    pattern = ''.join(tmp)
                    testing = df['Peptide amino acid sequence'].apply(motif_count, args=([pattern], 'cterm'))
                    top_ct = testing.iloc[:thresh].sum()
                    #bottom_ct = testing.iloc[thresh:-thresh].sum() + 1
                    bottom_ct = testing.sum()
                    top_n = thresh
                    #bot_n = len(testing)-thresh*2+1
                    bot_n = len(testing)
                    pval = stats.binom_test(top_ct, n=top_n, p=bottom_ct/bot_n, alternative='greater')
                    #output_list.append([pattern, i, j, aa1, aa2, top_ct, bottom_ct, thresh, len(testing)-2*thresh, pval])
                    output_list.append([pattern, i, j, aa1, aa2, top_ct, bottom_ct, thresh, len(testing), pval])

    # compile results
    mycols = ['motif', 'pos1', 'pos2', 'aa1', 'aa2', 'top ct', 'bot ct', 'top total', 'bot total', 'pvalue']
    result = pd.DataFrame(output_list, columns=mycols)
    result['log(OR)'] = np.log2((result['top ct'] / (result['top total']-result['top ct'])) / (result['bot ct'] / (result['bot total']-result['bot ct'])))
    result['adjusted pvalue'] = result['pvalue'] * len(result)
    result['qvalue'] = pvalue.bh_fdr(result['pvalue'])
    result.sort_values('pvalue', inplace=True)

    return result


def cterm_trimotif_enrichment(df, thresh):
    aa_vocab = degron_pred.vocab
    # get background probs for each aa
    aa_bg = {}
    #bot_n = len(df)-thresh*2+1
    bot_n = len(df)
    for aa in aa_vocab:
        tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([aa], 'cterm'))
        #tmp_prob = (tmp.iloc[thresh:-thresh].sum() + 1) / (bot_n*6)
        tmp_prob = (tmp.sum()) / (bot_n*6)
        aa_bg[aa] = tmp_prob

    output_list = []
    for aa1 in aa_vocab:
        for aa2 in aa_vocab:
            for aa3 in aa_vocab:
                pattern = ''.join([aa1, aa2, aa3])
                # figure out full count
                testing = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern], 'cterm'))
                top_ct = testing.iloc[:thresh].sum()
                top_n = thresh * 4
                # get baseline prob
                bottom_prob_1 = aa_bg[pattern[0]]
                tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern[:2]], 'cterm'))
                #bottom_prob_2_di = (tmp.iloc[thresh:-thresh].sum() + 1) / (bot_n*5)
                bottom_prob_2_di = (tmp.sum()) / (bot_n*5)
                bottom_prob_2 = bottom_prob_2_di / bottom_prob_1
                bottom_prob_2_aa = aa_bg[pattern[1]]
                tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern[1:3]], 'cterm'))
                #bottom_prob_3_di = (tmp.iloc[thresh:-thresh].sum() + 1) / (bot_n*5)
                bottom_prob_3_di = (tmp.sum()) / (bot_n*5)
                bottom_prob_3 = bottom_prob_3_di / bottom_prob_2_aa
                bottom_prob = bottom_prob_1 * bottom_prob_2 * bottom_prob_3

                # measure significance
                pval = stats.binom_test(top_ct, n=top_n, p=bottom_prob, alternative='greater')
                #output_list.append([pattern, top_ct, bottom_prob, thresh, len(testing)-2*thresh, pval])
                output_list.append([pattern, top_ct, bottom_prob, thresh, bot_n, pval])

    # compile results
    mycols = ['motif', 'top ct', 'background p', 'top total', 'bot total', 'pvalue']
    result = pd.DataFrame(output_list, columns=mycols)
    result['log(OR)'] = np.log2((result['top ct'] / (result['top total']-result['top ct'])) / (result['background p'] / (1 - result['background p'])))
    result['adjusted pvalue'] = result['pvalue'] * len(result)
    result.loc[result['adjusted pvalue']>1, 'adjusted pvalue'] = 1
    result['qvalue'] = pvalue.bh_fdr(result['pvalue'])
    result.sort_values('pvalue', inplace=True)

    return result


def cterm_quadmotif_enrichment(df, thresh):
    aa_vocab = degron_pred.vocab
    # get background probs for each aa
    aa_bg = {}
    #bot_n = len(df)-thresh*2+1
    bot_n = len(df)
    for aa in aa_vocab:
        tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([aa], 'cterm'))
        #tmp_prob = (tmp.iloc[thresh:-thresh].sum() + 1) / (bot_n*6)
        tmp_prob = (tmp.sum()) / (bot_n*6)
        aa_bg[aa] = tmp_prob

    output_list = []
    for aa1 in aa_vocab:
        for aa2 in aa_vocab:
            for aa3 in aa_vocab:
                for aa4 in aa_vocab:
                    pattern = ''.join([aa1, aa2, aa3, aa4])
                    # figure out full count
                    testing = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern], 'cterm'))
                    top_ct = testing.iloc[:thresh].sum()
                    top_n = thresh * 3
                    # get baseline prob
                    bottom_prob_1 = aa_bg[pattern[0]]
                    tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern[:2]], 'cterm'))
                    bottom_prob_2_di = (tmp.sum()) / (bot_n*5)
                    bottom_prob_2 = bottom_prob_2_di / bottom_prob_1
                    bottom_prob_2_aa = aa_bg[pattern[1]]
                    tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern[1:3]], 'cterm'))
                    bottom_prob_3_di = (tmp.sum()) / (bot_n*5)
                    bottom_prob_3 = bottom_prob_3_di / bottom_prob_2_aa
                    bottom_prob_3_aa = aa_bg[pattern[2]]
                    tmp = df['Peptide amino acid sequence'].apply(motif_count_nopos, args=([pattern[2:4]], 'cterm'))
                    bottom_prob_4_di = (tmp.sum()) / (bot_n*5)
                    bottom_prob_4 = bottom_prob_4_di / bottom_prob_3_aa

                    # calc background p
                    bottom_prob = bottom_prob_1 * bottom_prob_2 * bottom_prob_3 * bottom_prob_4

                    # measure significance
                    pval = stats.binom_test(top_ct, n=top_n, p=bottom_prob, alternative='greater')
                    #output_list.append([pattern, top_ct, bottom_prob, thresh, len(testing)-2*thresh, pval])
                    output_list.append([pattern, top_ct, bottom_prob, thresh, bot_n, pval])

    # compile results
    mycols = ['motif', 'top ct', 'background p', 'top total', 'bot total', 'pvalue']
    result = pd.DataFrame(output_list, columns=mycols)
    result['log(OR)'] = np.log2((result['top ct'] / (result['top total']-result['top ct'])) / (result['background p'] / (1 - result['background p'])))
    result['adjusted pvalue'] = result['pvalue'] * len(result)
    result.loc[result['adjusted pvalue']>1, 'adjusted pvalue'] = 1
    result['qvalue'] = pvalue.bh_fdr(result['pvalue'])
    result.sort_values('pvalue', inplace=True)

    return result


def main(opts):
    # read data
    df = pd.read_csv(opts['input'], sep='\t').rename(columns={'regulatory potential': 'degron potential'})
    df.sort_values('degron potential', ascending=False, inplace=True)

    # mark sequences with class lables
    df['y'] = 0
    if opts['type'] == 'cterm':
        df.loc[df['Modal Bin']<2.5, 'y'] = 1
    else:
        df.loc[df['Modal Bin']<3.5, 'y'] = 1

    thresh = calc_top_degron_threshold(df)
    if opts['type'] == 'nterm':
        # analyze enrich for n-terminal motif
        di_result = nterm_dimotif_enrichment(df.copy(), thresh)
        tri_result = nterm_trimotif_enrichment(df.copy(), thresh)
    else:
        # analyze enrichment for c-terminal motifs
        di_result = cterm_dimotif_enrichment(df.copy(), thresh)
        tri_result = cterm_trimotif_enrichment(df.copy(), thresh)
        #quad_result = cterm_quadmotif_enrichment(df.copy(), thresh)

    # save results
    di_path = os.path.join(opts['output'], 'dimotif.txt')
    di_result.to_csv(di_path, sep='\t', index=False)
    tri_path = os.path.join(opts['output'], 'trimotif.txt')
    tri_result.to_csv(tri_path, sep='\t', index=False)
    #quad_path = os.path.join(opts['output'], 'quadmotif.txt')
    #quad_result.to_csv(quad_path, sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
