import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
import deepDegron.statistical_test as test


def test_mn1_main():
    # make output directory if it doesn't exist
    outdir = os.path.join(file_dir, 'output')
    if not os.path.exists(outdir): os.mkdir(outdir)

    cterm_model1 = os.path.join(file_dir, "../models/cterm/neural_network_pos_specific.pickle")
    cterm_model2 = os.path.join(file_dir, "../models/cterm/neural_network_bag_of_amino_acids.pickle")

    # run ctnnb1 example
    opts = {
        'input': os.path.join(file_dir, 'data/MN1_example.txt'),
        'degrons': None,
        'sites': None, #os.path.join(file_dir, 'data/ctnnb1_phosphodegron_sites.txt'),
        'flank': 0,
        'num_sim': 10000,
        'processes': 0,
        'cterm_models': '{0},{1}'.format(cterm_model1, cterm_model2),
        'ensembl_release': 75,
        'nterm_models': None,
        'truncation': None,
        'output': os.path.join(file_dir, 'output/MN1_results.txt')
    }
    test.main(opts)
