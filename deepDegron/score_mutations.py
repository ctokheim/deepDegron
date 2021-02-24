"""
File: score_mutations.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Script to score degron potential for mutations
"""
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
os.environ['OMP_NUM_THREADS'] = "1"

import argparse
import pandas as pd
import numpy as np
import varcode
from pyensembl import ensembl_grch37, EnsemblRelease
import deepDegron.sequence_context as sc
from varcode.effects import *
import collections
import deepDegron.utils as utils
import deepDegron.degron_pred as degron_pred
from deepDegron.variants import *
import deepDegron.simulation as simulation
import deepDegron.pvalue as pvalue
import csv

# logging import
import logging
logger = logging.getLogger(__name__)  # module logger
logging.getLogger('varcode').setLevel(logging.ERROR)
logging.getLogger('tensorflow').setLevel(logging.ERROR)
logging.getLogger('Bio').setLevel(logging.ERROR)
import warnings
warnings.filterwarnings("ignore")
# multiprocess
from multiprocessing import Pool

def parse_arguments():
    info = 'Script to score mutations for degron potential'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutation data in MAF format')
    parser.add_argument('-c', '--cterm-models',
                        type=str, default=None,
                        help='Path to saved cterminal degron models')
    parser.add_argument('-n', '--nterm-models',
                        type=str, default=None,
                        help='Path to saved nterminal degron models')
    parser.add_argument('-e', '--ensembl-release',
                        type=int, default=75,
                        help='Ensembl release version number for gene annotations in varcode (default: 75)')
    parser.add_argument('-r', '--raw',
                        action='store_true', default=False,
                        help='Use raw deepDegron score')
    parser.add_argument('-p', '--processes',
                        type=int, default=0,
                        help='Number of computer processes to parallelize calculations')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Result file')
    args = parser.parse_args()
    return vars(args)

def check_analysis_type(options):
    """Determines what type of analysis to do based on the command line options."""
    if options['cterm_models']:
        analysis_type = 'cterminus'
    else:
        analysis_type = 'nterminus'
    return analysis_type

def singleprocess_score(info):
    """Unpacks the multiprocess input"""
    options, mychr = info
    analysis_type = check_analysis_type(options)
    return analyze(options, chrom=mychr, analysis=analysis_type, release=options['ensembl_release'])

def multiprocess_score(opts):
    """Handles parallelization of permutations by splitting work
    by chromosome.
    """
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    chroms = list(map(str, range(1, 23))) + ['X', 'Y']
    result_list = []
    if multiprocess_flag:
        for i in range(0, len(chroms), num_processes):
            pool = Pool(processes=num_processes)
            tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
            info_repeat = ((opts, chroms[tmp_ix])
                            for tmp_ix in range(i, i+tmp_num_proc))
            process_results = pool.imap(singleprocess_score, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                for chrom_result in process_results:
                    if chrom_result:
                        result_list += chrom_result
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
    else:
        analysis_type = check_analysis_type(opts)
        result_list += analyze(opts, analysis=analysis_type, release=opts['ensembl_release'])

    return result_list

def score(variant_list,
          tx, clf1, clf2,
          model='cterm',
          nuc_context=1.5):
    """Simulate the effect of mutations on terminal degrons. Handles both c-terminal and and n-terminal degrons.

    """
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    trinuc_context = [sc.get_chasm_context(nc) for nc in trinuc_context]
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides
    var_sub_no_nmd = utils.filter_nmd_subs(var_sub, tx)

    # filter for indel variants
    var_indel = [x for x in variant_list
                 if x.__class__ in utils.indels
                 and x.variant.is_indel]
    dna_change_indel = [[v.variant.contig, v.variant.start, v.variant.ref, v.variant.alt]
                        for v in var_indel]
    var_indel = utils.nmd(var_indel, tx, drop=True)

    # return if no substitution
    if not var_sub and not var_indel:
        return []

    # figure out the affect on cterminal degrons
    myvars, delta_prob = degron_pred.delta_prob(var_sub_no_nmd+var_indel, tx, clf1, clf2, model=model, is_sum=False)

    # skip if no terminal variants
    if not len(delta_prob):
        return []

    # fill in info about each variant
    var_info = []
    for ix, var in enumerate(myvars):
        tmp = [
            var.gene_name, var.transcript_id,
            var.variant.contig, var.variant.start, var.variant.end,
            var.short_description,
            delta_prob.iloc[ix]
        ]
        var_info.append(tmp)

    return var_info

def score_raw(variant_list,
          tx, clf1, clf2,
          model='cterm',
          nuc_context=1.5):
    """Simulate the effect of mutations on terminal degrons. Handles both c-terminal and and n-terminal degrons.

    """
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    trinuc_context = [sc.get_chasm_context(nc) for nc in trinuc_context]
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides
    var_sub_no_nmd = utils.filter_nmd_subs(var_sub, tx)

    # filter for indel variants
    var_indel = [x for x in variant_list
                 if x.__class__ in utils.indels
                 and x.variant.is_indel]
    dna_change_indel = [[v.variant.contig, v.variant.start, v.variant.ref, v.variant.alt]
                        for v in var_indel]
    # skip dropping based on nmd
    #var_indel = utils.nmd(var_indel, tx, drop=True)

    # return if no substitution
    if not var_sub and not var_indel:
        return []

    # figure out the affect on cterminal degrons
    myvars, delta_prob, prob = degron_pred.delta_prob_raw(var_sub_no_nmd+var_indel, tx, clf1, clf2, model=model, is_sum=False)

    # skip if no terminal variants
    if not len(delta_prob):
        return []

    # fill in info about each variant
    var_info = []
    for ix, var in enumerate(myvars):
        tmp = [
            var.gene_name, var.transcript_id,
            var.variant.contig, var.variant.start, var.variant.end,
            var.short_description,
            delta_prob.iloc[ix], prob.iloc[ix]
        ]
        var_info.append(tmp)

    return var_info

def analyze(opts, chrom=None, release=75, analysis='cterminus'):
    # read in data
    variant_dict = read_maf(opts['input'], chrom, release)

    # read the c-terminus classifier
    if analysis == 'cterminus':
        clf1_path, clf2_path = opts['cterm_models'].split(',')
        clf1 = degron_pred.load_classifier(clf1_path)
        clf2 = degron_pred.load_classifier(clf2_path)
    # read the n-terminus classifier
    if analysis == 'nterminus':
        clf1_path, clf2_path = opts['nterm_models'].split(',')
        clf1 = degron_pred.load_classifier(clf1_path)
        clf2 = degron_pred.load_classifier(clf2_path)

    # iterate over each gene in the MAF file
    if not opts['raw']:
        output_list = [['Hugo_Symbol', 'Transcript_ID', 'Chromosome',
                        'Start_Position', 'End_Position', 'HGVSp_Short',
                        'delta degron potential']]
    else:
        output_list = [['Hugo_Symbol', 'Transcript_ID', 'Chromosome',
                        'Start_Position', 'End_Position', 'HGVSp_Short',
                        'delta deepDegron score', 'deepDegron score']]
    for gene in variant_dict:
        # variants for a specific gene
        ensembl_tx_name = variant_dict[gene]['transcript_name']
        tx = variant_dict[gene]['transcript']
        variant_list = variant_dict[gene]['variants']

        # skip non-protein coding genes
        if tx.biotype != 'protein_coding' or not tx.complete:
            continue

        # calculate the significance
        if not opts['raw']:
            if analysis == 'cterminus':
                results = score(variant_list, tx, clf1, clf2, model='cterm')
            elif analysis == 'nterminus':
                results = score(variant_list, tx, clf1, clf2, model='nterm')
        else:
            if analysis == 'cterminus':
                results = score_raw(variant_list, tx, clf1, clf2, model='cterm')
            elif analysis == 'nterminus':
                results = score_raw(variant_list, tx, clf1, clf2, model='nterm')

        # append results
        if results:
            output_list += list(results)

    return output_list

def main(opts):
    result = multiprocess_score(opts)
    with open(opts['output'], 'w') as whandle:
        mywriter = csv.writer(whandle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(result)

def cli_main():
    # run main with CLI options
    opts = parse_arguments()
    main(opts)


if __name__ == '__main__':
    cli_main()
