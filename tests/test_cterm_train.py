import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['KERAS_BACKEND'] = 'tensorflow'
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
import deepDegron.train_cterm_nn as train


def test_train_cterm():
    # make output directory if it doesn't exist
    outdir = os.path.join(file_dir, 'output')
    if not os.path.exists(outdir): os.mkdir(outdir)

    # train c-terminal model
    opts = {
        'input': os.path.join(file_dir, 'data/gps_cterminal_degron_screen_1k.txt'),
        'validation_perf': os.path.join(file_dir, 'output/cterm_validation_performance.txt'),
        'test_pred': os.path.join(file_dir, 'output/cterm_test_set_prediction.txt'),
        'bag_of_aa': os.path.join(file_dir, 'output/cterm_neural_network_bag_of_amino_acids.pickle'),
        'sequence_specific': os.path.join(file_dir, 'output/cterm_neural_network_pos_specific.pickle'),
        'output': os.path.join(file_dir, 'output/cterminal_degron_predictions.txt')
    }
    try:
        train.main(opts)
    except:
        import IPython ; IPython.embed() ; raise

test_train_cterm()
