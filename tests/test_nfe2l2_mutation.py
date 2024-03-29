import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
import deepDegron.statistical_test as test


def test_nfe2l2_main():
    # make output directory if it doesn't exist
    outdir = os.path.join(file_dir, 'output')
    if not os.path.exists(outdir): os.mkdir(outdir)

    # run nfe2l2 example
    opts = {
        'input': os.path.join(file_dir, 'data/NFE2L2_mutations_lusc.txt'),
        'degrons': os.path.join(file_dir, 'data/nfe2l2_degron_sites.txt'),
        'sites': None,
        'flank': 0,
        'num_sim': 100,
        'processes': 0,
        'cterm_models': None,
        'nterm_models': None,
        'ensembl_release': 75,
        'output': os.path.join(file_dir, 'output/nfe2l2_degron_results.txt')
    }
    test.main(opts)
