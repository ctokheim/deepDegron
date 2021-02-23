import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
import deepDegron.score_mutations as score


def test_score_main():
    # make output directory if it doesn't exist
    outdir = os.path.join(file_dir, 'output')
    if not os.path.exists(outdir): os.mkdir(outdir)

    cterm_model1 = os.path.join(file_dir, "../models/cterm/neural_network_pos_specific.pickle")
    cterm_model2 = os.path.join(file_dir, "../models/cterm/neural_network_bag_of_amino_acids.pickle")

    # score gata3 mutations
    opts = {
        'input': os.path.join(file_dir, 'data/gata3_mutations.txt'),
        'raw': True,
        'processes': 0,
        'cterm_models': '{0},{1}'.format(cterm_model1, cterm_model2),
        'nterm_models': None,
        'ensembl_release': 75,
        'output': os.path.join(file_dir, 'output/gata3_delta_degron_scores.txt')
    }
    score.main(opts)
