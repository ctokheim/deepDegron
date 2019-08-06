import pandas as pd
import numpy as np
import scipy.stats as stats
import pickle
import degron_pred

import argparse


def parse_arguments():
    info = 'Score the saturation mutagenesis results'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Saturation mutagenesis file')
    parser.add_argument('-n', '--neural-network',
                        type=str, required=True,
                        help='Two neural networks (position specific and bag of amino acids separated by comma)')
    parser.add_argument('-t', '--type',
                        type=str, required=True,
                        help='Type of model (cterm or nterm)')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Prediction results of saturation mutagenesis')
    args = parser.parse_args()
    return vars(args)


def compute_cterm_feature_matrix(sequences, split, di=False):
    """Compute the feature matrix for c-terminal degrons"""
    if 0 < split < 23:
        X = degron_pred.binned_bag_of_words(sequences.str[-split:],
                                            int(split), n=int(split),
                                            dinuc=di)
        X2 = degron_pred.binned_bag_of_words(sequences.str[:-split],
                                             1, n=23-int(split), dinuc=False)
        X = np.hstack([X2.toarray(), X])
    elif split == 0:
        X = degron_pred.binned_bag_of_words(sequences,
                                            int(split), n=int(split),
                                            dinuc=False)
    elif split == 23:
        X = degron_pred.binned_bag_of_words(sequences,
                                            int(split), n=int(split),
                                            dinuc=di)
    return X


def compute_nterm_feature_matrix(sequences, split, di=False):
    """Compute the feature matrix for n-terminal degrons"""
    if 0 < split < 23:
        X = degron_pred.binned_bag_of_words(sequences.str[:split],
                                            int(split), n=int(split),
                                            dinuc=di, cterm=False)
        X2 = degron_pred.binned_bag_of_words(sequences.str[split:],
                                             1, n=23-int(split),
                                             dinuc=False, cterm=False)
        X = np.hstack([X2.toarray(), X])
    elif split == 0:
        X = degron_pred.binned_bag_of_words(sequences,
                                            int(split), n=int(split),
                                            dinuc=False, cterm=False)
    elif split == 23:
        X = degron_pred.binned_bag_of_words(sequences,
                                            int(split), n=int(split),
                                            dinuc=di, cterm=False)
    return X


def main(opts):
    # read in classifiers
    pos_path, bag_path = opts['neural_network'].split(',')
    with open(pos_path, 'rb') as handle:
        seq_clf = pickle.load(handle)
    with open(bag_path, 'rb') as handle:
        bag_clf = pickle.load(handle)

    # read in sat mut data file
    sat_mut = pd.read_csv(opts['input'], sep='\t')

    seq_col = 'Peptide amino acid sequence'
    if opts['type'] == 'cterm':
        X_seq = compute_cterm_feature_matrix(sat_mut[seq_col], 6, di=True)
        X_bag = compute_cterm_feature_matrix(sat_mut[seq_col], 0, di=False)
    else:
        X_seq = compute_nterm_feature_matrix(sat_mut[seq_col], 6, di=True)
        X_bag = compute_nterm_feature_matrix(sat_mut[seq_col], 0, di=False)

    # score the predictions
    seq_pred = seq_clf.predict_proba(X_seq)[:,0]
    bag_pred = bag_clf.predict_proba(X_bag)[:,0]

    # add info to DF
    sat_mut['sequence position specific'] = seq_pred
    sat_mut['bag of words'] = bag_pred
    sat_mut['degron potential'] = sat_mut['sequence position specific'] - sat_mut['bag of words']

    # save results
    sat_mut.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
