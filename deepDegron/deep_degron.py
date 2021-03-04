"""
File: deep_degron.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Command line parsing script for package
"""

# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
os.environ['OMP_NUM_THREADS'] = "1"
os.environ['KMP_WARNINGS'] = '0'
import argparse
import warnings
warnings.filterwarnings("ignore")


def parse_arguments():
    info = 'deepDegron command line tool to predict degrons'
    parent_parser = argparse.ArgumentParser(description=info)

    # add subparsers
    subparsers = parent_parser.add_subparsers(title='Sub-commands', dest='kind')
    help_info = 'Score mutations for predicted impact on degrons.'
    parser_score = subparsers.add_parser('score',
                                      help=help_info,
                                      description=help_info + 'DeepDegron can score a protein sequeqnce for potential to contain a degron.'
                                      'The score sub-command assesses whether a mutation changes deepDegron\'s prediction of a degron '
                                      '(i.e. delta degron potential).')
    help_info = 'Evaluates whether a gene contains a signficant enrichment of mutations prediction to impact a degron.'
    parser_test = subparsers.add_parser('test',
                                        help=help_info,
                                        description='The "test" sub-command performs a statistical test '
                                        'to evaluate whether mutations preferentially disrupt a predicted degron in a gene.')
    help_info = 'Train a neural network model for c-terminal degrons.'
    parser_cterm = subparsers.add_parser('cterm_train',
                                         help=help_info,
                                         description=help_info)
    help_info = 'Train a neural network model for n-terminal degrons.'
    parser_nterm = subparsers.add_parser('nterm_train',
                                         help=help_info,
                                         description=help_info)

    ###################
    # specify command line options for each sub-parser
    ###################
    # score sub-parser
    parser_score.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutation data in MAF format')
    parser_score.add_argument('-c', '--cterm-models',
                        type=str, default=None,
                        help='Path to saved cterminal degron models')
    parser_score.add_argument('-n', '--nterm-models',
                        type=str, default=None,
                        help='Path to saved nterminal degron models')
    parser_score.add_argument('-e', '--ensembl-release',
                        type=int, default=75,
                        help='Ensembl release version number for gene annotations in varcode (default: 75)')
    parser_score.add_argument('-r', '--raw',
                        action='store_true', default=False,
                        help='Use raw deepDegron score')
    parser_score.add_argument('-p', '--processes',
                        type=int, default=0,
                        help='Number of computer processes to parallelize calculations')
    parser_score.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Result file')

    # test sub-parser
    parser_test.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutation data in MAF format')
    parser_test.add_argument('-d', '--degrons',
                        type=str, default=None,
                        help='Degron locations')
    parser_test.add_argument('-t', '--truncation',
                        action='store_true', default=False,
                        help='Analyze clustering of truncating mutations')
    parser_test.add_argument('-s', '--sites',
                        type=str, default=None,
                        help='Sites of interest')
    parser_test.add_argument('-c', '--cterm-models',
                        type=str, default=None,
                        help='Path to saved cterminal degron models')
    parser_test.add_argument('-n', '--nterm-models',
                        type=str, default=None,
                        help='Path to saved nterminal degron models')
    parser_test.add_argument('-f', '--flank',
                        type=int, default=0,
                        help='Number of flanking residues to consider')
    parser_test.add_argument('-ns', '--num-sim',
                        type=int, default=10000,
                        help='Number of simulations in statistical test')
    parser_test.add_argument('-p', '--processes',
                        type=int, default=0,
                        help='Number of computer processes to parallelize calcuations')
    parser_test.add_argument('-e', '--ensembl-release',
                        type=int, default=75,
                        help='Ensembl release version number for gene annotations in varcode (default: 75)')
    parser_test.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Result file')

    # cterm_train sub-parser
    parser_cterm.add_argument('-i', '--input',
                        type=str, required=True,
                        help='training data')
    parser_cterm.add_argument('-v', '--validation-perf',
                        type=str, required=True,
                        help='Validation performance')
    parser_cterm.add_argument('-t', '--test-pred',
                        type=str, required=True,
                        help='File to save predictions for the test set')
    parser_cterm.add_argument('-b', '--bag-of-aa',
                        type=str, required=True,
                        help='bag of amino acids trained model')
    parser_cterm.add_argument('-s', '--sequence-specific',
                        type=str, required=True,
                        help='sequence specific trained model')
    parser_cterm.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Predictions for every sequence')

    # nterm_train sub-parser
    parser_nterm.add_argument('-i', '--input',
                        type=str, required=True,
                        help='training data')
    parser_nterm.add_argument('-v', '--validation-perf',
                        type=str, required=True,
                        help='Validation performance')
    parser_nterm.add_argument('-t', '--test-pred',
                        type=str, required=True,
                        help='File to save predictions for the test set')
    parser_nterm.add_argument('-b', '--bag-of-aa',
                        type=str, required=True,
                        help='bag of amino acids trained model')
    parser_nterm.add_argument('-s', '--sequence-specific',
                        type=str, required=True,
                        help='sequence specific trained model')
    parser_nterm.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Predictions for every sequence')

    args = parent_parser.parse_args()
    return vars(args)


def main(opts):
    if opts['kind'] == 'score':
        import deepDegron.score_mutations as score
        score.main(opts)
    elif opts['kind'] == 'test':
        import deepDegron.statistical_test as test
        test.main(opts)
    elif opts['kind'] == 'cterm_train':
        import deepDegron.train_cterm_nn as ctrain
        ctrain.main(opts)
    elif opts['kind'] == 'nterm_train':
        import deepDegron.train_nterm_nn as ntrain
        ntrain.main(opts)


def cli_main():
    # run main with CLI options
    opts = parse_arguments()
    main(opts)


if __name__ == '__main__':
    cli_main()

