"""
File: statistical_test.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Build null distribution
"""
import argparse
import pandas as pd
import numpy as np
import varcode
from pyensembl import ensembl_grch37, EnsemblRelease
import sequence_context as sc
from varcode.effects import *
import collections
import utils
import degron_pred
from variants import *
import simulation
import pvalue

# logging import
import logging
logger = logging.getLogger(__name__)  # module logger
# multiprocess
from multiprocessing import Pool

def parse_arguments():
    info = 'Statistical test'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutation data in MAF format')
    parser.add_argument('-d', '--degrons',
                        type=str, default=None,
                        help='Degron locations')
    parser.add_argument('-s', '--sites',
                        type=str, default=None,
                        help='Sites of interest')
    parser.add_argument('-c', '--cterm-models',
                        type=str, default=None,
                        help='Path to saved cterminal degron models')
    parser.add_argument('-f', '--flank',
                        type=int, default=0,
                        help='Number of flanking residues to consider')
    parser.add_argument('-p', '--processes',
                        type=int, default=0,
                        help='Number of processes')
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
        #analysis_type = 'degrons' if opts['degrons'] else 'lysine'
        analysis_type = 'cterminus'
        result_list += analyze(opts, analysis=analysis_type)

    return result_list


def singleprocess_permutation(info):
    """Unpacks the multiprocess input"""
    options, mychr = info
    if options['degrons']:
        analysis_type = 'degrons'
    elif options['sites']:
        analysis_type = 'sites'
    elif options['cterm_models']:
        analysis_type = 'cterminus'
    else:
        analysis_type = 'lysine'
    return analyze(opts, chrom=mychr, analysis=analysis_type)


def analyze(opts, chrom=None, analysis='degrons'):
    # read in data
    variant_dict = read_maf(opts['input'], chrom)

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
        #clf1 = degron_pred.load_classifier('data/logistic_regression_pos_specific_dinuc_v2.pickle')
        #clf2 = degron_pred.load_classifier('data/logistic_regression_bag_of_words_v2.pickle')

    # iterate over each gene in the MAF file
    output_list = []
    for gene in variant_dict:
        # variants for a specific gene
        ensembl_tx_name = variant_dict[gene]['transcript_name']
        tx = variant_dict[gene]['transcript']
        variant_list = variant_dict[gene]['variants']

        # only consider genes with degrons
        if analysis == 'degrons' and gene not in degron_intvls:
            continue

        # only consider genes with ub sites
        if analysis == 'sites' and ensembl_tx_name not in ub_intvls:
            continue

        # skip non-protein coding genes
        if tx.biotype != 'protein_coding' or not tx.complete:
            continue

        # calculate the significance
        if analysis == 'degrons':
            results = simulation.degron(variant_list, tx, degron_intvls[gene])
        elif analysis == 'cterminus':
            results = simulation.cterm_degron(variant_list, tx, clf1, clf2)
        elif analysis == 'sites':
            results = simulation.site(variant_list, tx, ub_intvls[ensembl_tx_name])
        else:
            results = simulation.lysine_mutations(variant_list, tx)

        # append results
        if results:
            output_list.append([gene]+list(results))

    return output_list
    """
    mycols = ['gene', 'num_degrons_mut', 'pval']
    output_df = pd.DataFrame(output_list, columns=mycols)

    import IPython ; IPython.embed()
    raise

    mycols = ['gene', 'new_lys', 'new_lys_pval', 'lost_lys', 'lost_lys_pval']
    output_df = pd.DataFrame(output_list, columns=mycols)

    import IPython ; IPython.embed()
    raise
    """

    """
    ####################
    # Look into indels
    ####################
    var_indel = [x for x in variant_list
                 if x.__class__ in utils.indels
                 and x.variant.is_indel]
    dna_change_indel = [[v.variant.contig, v.variant.start, v.variant.ref, v.variant.alt]
                        for v in var_indel]

    # there is no sequence context for indels
    sc_indel = sc.SequenceContext(gata3_tx, 0)

    tmp_input = [('None', len(var_indel))]
    random_indel = sc_indel.random_pos(tmp_input, 10000)

    sim_delta_deg_list = []
    tmp_mut_pos = np.hstack(abs_pos for ctxt, cds_pos, abs_pos in random_indel)
    for sim_pos in tmp_mut_pos:
        # get variant effects for simulated mutations
        sim_variant_indels = get_mutation_info(sim_pos, gata3_tx, dna_change_indel)

        # predict whether mutation will be sensitive to NMD by rule of thumb
        sim_nmd_sensitivity = utils.nmd(sim_variant_indels, gata3_tx)

        # prepare df
        cterm_seq = [v.mutant_protein_sequence[-23:]+'*'
                     if v.mutant_protein_sequence
                     else 'A'*23+'*'
                     for v in sim_variant_indels]
        sim_result_df = pd.DataFrame({'seq': cterm_seq, 'nmd': sim_nmd_sensitivity})

        # create feature matrix
        X = degron_pred.binned_bag_of_words(sim_result_df['seq'], 23, dinuc=True)
        X2 = degron_pred.binned_bag_of_words(sim_result_df['seq'], 1)
        sim_result_df['prob'] = clf.predict_proba(X)[:, 1]
        sim_result_df['prob2'] = clf2.predict_proba(X2)[:, 1]
        sim_result_df['delta prob'] = sim_result_df['prob'] - sim_result_df['prob2']

        # now compute the test statistic of interest
        #sim_delta_deg = (sim_result_df.loc[sim_result_df['nmd']==0, 'prob']-0.7859).sum()
        sim_delta_deg = (sim_result_df.loc[sim_result_df['nmd']==0, 'delta prob']-0.32).sum()
        sim_delta_deg_list.append(sim_delta_deg)


    # look at score distribution for actually observed frameshifts
    nmd_sensitivity = utils.nmd(var_indel, gata3_tx)
    cterm_seq = [v.mutant_protein_sequence[-23:]+'*'
                 if v.mutant_protein_sequence
                 else 'A'*23+'*'
                 for v in var_indel]
    result_df = pd.DataFrame({'seq': cterm_seq, 'nmd': nmd_sensitivity})

    # create feature matrix
    X = degron_pred.binned_bag_of_words(result_df['seq'], 23, dinuc=True)
    X2 = degron_pred.binned_bag_of_words(result_df['seq'], 1)
    result_df['prob'] = clf.predict_proba(X)[:, 1]
    result_df['prob2'] = clf2.predict_proba(X2)[:, 1]
    result_df['delta prob'] = result_df['prob'] - result_df['prob2']

    # now compute the test statistic of interest
    #delta_deg = (result_df.loc[result_df['nmd']==0, 'prob']-0.7859).sum()
    delta_deg = (result_df.loc[result_df['nmd']==0, 'delta prob']-0.32).sum()

    import IPython ; IPython.embed()
    raise


    #sim_variant_indels = []
    #for i in range(len(sim1)):
        #contig, start, ref, alt = dna_change_indel[i]
        #new_pos = sim1[i]

        # figure out the new reference seq
        # for deletions
        #len_del = len(ref)
        #if len_del:
            #tmp_offset_pos = gata3_tx.spliced_offset(int(new_pos)) - start_offset
            #new_ref = gata3_tx.coding_sequence[tmp_offset_pos:tmp_offset_pos+len_del]

        #new_var_info = [contig, new_pos, new_ref, alt]
        #sim_variant_indels.append(read_variant(new_var_info, gata3_tx))

    #Variant(contig=7, start=140453136, ref="A", alt="T", ensembl=ensembl_grch37)

    import IPython ; IPython.embed()
    """


def main(opts):
    # do the analysis
    results = multiprocess_permutation(opts)

    # format the output
    if opts['degrons']:
        output_df = utils.process_degron_results(results)
    elif opts['sites']:
        output_df = utils.process_ub_results(results)
    else:
        output_df = utils.process_cterm_degron_results(results)

    # save the results
    output_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
