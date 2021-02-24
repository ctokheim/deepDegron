"""
File: statistical_test.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Perform a statistical test of whether there is enrichment of mutations at degrons
"""
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
os.environ['OMP_NUM_THREADS'] = "1"
os.environ['KMP_WARNINGS'] = '0'

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
    info = 'Performs a statistical test for the enrichment of mutations at degrons'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutation data in MAF format')
    parser.add_argument('-d', '--degrons',
                        type=str, default=None,
                        help='Degron locations')
    parser.add_argument('-t', '--truncation',
                        action='store_true', default=False,
                        help='Analyze clustering of truncating mutations')
    parser.add_argument('-s', '--sites',
                        type=str, default=None,
                        help='Sites of interest')
    parser.add_argument('-c', '--cterm-models',
                        type=str, default=None,
                        help='Path to saved cterminal degron models')
    parser.add_argument('-n', '--nterm-models',
                        type=str, default=None,
                        help='Path to saved nterminal degron models')
    parser.add_argument('-f', '--flank',
                        type=int, default=0,
                        help='Number of flanking residues to consider')
    parser.add_argument('-ns', '--num-sim',
                        type=int, default=10000,
                        help='Number of simulations in statistical test')
    parser.add_argument('-p', '--processes',
                        type=int, default=0,
                        help='Number of computer processes to parallelize calcuations')
    parser.add_argument('-e', '--ensembl-release',
                        type=int, default=75,
                        help='Ensembl release version number for gene annotations in varcode (default: 75)')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Result file')
    args = parser.parse_args()
    return vars(args)


def multiprocess_permutation(opts):
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
            process_results = pool.imap(singleprocess_permutation, info_repeat)
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
        result_list += analyze(opts, analysis=analysis_type)

    return result_list


def check_analysis_type(options):
    """Determines what type of analysis to do based on the command line options."""
    if options['degrons']:
        analysis_type = 'degrons'
    elif options['sites']:
        analysis_type = 'sites'
    elif options['cterm_models']:
        analysis_type = 'cterminus'
    elif options['nterm_models']:
        analysis_type = 'nterminus'
    elif options['truncation']:
        analysis_type = 'truncation'
    else:
        analysis_type = 'sites'
    return analysis_type


def singleprocess_permutation(info):
    """Unpacks the multiprocess input"""
    options, mychr = info
    analysis_type = check_analysis_type(options)
    return analyze(options, chrom=mychr, analysis=analysis_type)


def analyze(opts, chrom=None, analysis='degrons'):
    # read in data
    variant_dict = read_maf(opts['input'], chrom)#, release=opts['ensembl_release'])

    # read in the degron data
    if analysis == 'degrons':
        degron_intvls = utils.read_degron_intervals(opts['degrons'])
    if analysis == 'sites':
        ub_intvls = utils.read_sites(opts['sites'], opts['flank'])
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
    output_list = []
    for gene in variant_dict:
        # variants for a specific gene
        ensembl_tx_name = variant_dict[gene]['transcript_name']
        tx = variant_dict[gene]['transcript']
        variant_list = variant_dict[gene]['variants']
        num_vars = len(variant_list)

        # compute the average position of the variants
        try:
            aa_len = len(tx.coding_sequence)/3
            avg_pos = np.mean([v.aa_mutation_start_offset/aa_len
                               for v in variant_list
                               if v.aa_mutation_start_offset is not None])
        except:
            avg_pos = np.nan

        ## Skip genes not relevant for a particular analysis ##
        # only consider genes with degrons
        if analysis == 'degrons' and gene not in degron_intvls:
            continue
        # only consider genes with ub sites
        if analysis == 'sites' and ensembl_tx_name not in ub_intvls:
            continue
        # skip non-protein coding genes
        if tx.biotype != 'protein_coding' or not tx.complete:
            continue
        # skip genes lacking NMD sensitive region
        if analysis == 'truncation':
            normal_prot_len = len(tx.protein_sequence)
            nmd_insens = utils.get_nmd_insensitive_len(tx)
            frac_nmd_insens = float(nmd_insens) / normal_prot_len
            if frac_nmd_insens >= 1:
                continue

        # calculate the significance
        if analysis == 'degrons':
            #results = simulation.degron(variant_list, tx, degron_intvls[gene], num_simulations=opts['num_sim'])
            results = simulation.degron_with_indel(variant_list, tx, degron_intvls[gene], num_simulations=opts['num_sim'])
        elif analysis == 'cterminus':
            results = simulation.terminal_degron(variant_list, tx, clf1, clf2, model='cterm', num_simulations=opts['num_sim'])
        elif analysis == 'nterminus':
            results = simulation.terminal_degron(variant_list, tx, clf1, clf2, model='nterm', num_simulations=opts['num_sim'])
        elif analysis == 'sites':
            results = simulation.site(variant_list, tx, ub_intvls[ensembl_tx_name], num_simulations=opts['num_sim'])
        elif analysis == 'truncation':
            results = simulation.clustered_truncation(variant_list, tx, num_simulations=opts['num_sim'])

        # append results
        if results:
            output_list.append([gene]+list(results)+[num_vars, avg_pos])

    return output_list


def main(opts):
    # do the analysis
    results = multiprocess_permutation(opts)

    # format the output
    if opts['degrons']:
        output_df = utils.process_degron_results(results)
    elif opts['sites']:
        output_df = utils.process_ub_results(results)
    elif opts['cterm_models'] or opts['nterm_models']:
        output_df = utils.process_terminal_degron_results(results)
    elif opts['truncation']:
        output_df = utils.process_trunc_results(results)

    # save the results
    output_df.to_csv(opts['output'], sep='\t', index=False)


def cli_main():
    # run main with CLI options
    opts = parse_arguments()
    main(opts)


if __name__ == '__main__':
    cli_main()
