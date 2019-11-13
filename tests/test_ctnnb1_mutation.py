import deepDegron.statistical_test as test
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))


def test_ctnnb1_main():
    # make output directory if it doesn't exist
    outdir = os.path.join(file_dir, 'output')
    if not os.path.exists(outdir): os.mkdir(outdir)

    # run ctnnb1 example
    opts = {
        'input': os.path.join(file_dir, 'data/CTNNB1_mutations_ucec.txt'),
        'degrons': None,
        'sites': os.path.join(file_dir, 'data/ctnnb1_phosphodegron_sites.txt'),
        'flank': 3,
        'num_sim': 100,
        'processes': 0,
        'cterm_models': None,
        'nterm_models': None,
        'output': os.path.join(file_dir, 'output/ctnnb1_phospho_results.txt')
    }
    test.main(opts)