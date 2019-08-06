"""
File: train_nterm_nn.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Train n-terminal degron model
"""
import pandas as pd
import argparse
import degron_pred
from sklearn.model_selection import train_test_split
import sklearn.metrics as metrics
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
import pickle
import numpy as np
from sklearn.utils import class_weight


def parse_arguments():
    info = 'Train a c-terminal degron model'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='training data')
    parser.add_argument('-v', '--validation-perf',
                        type=str, required=True,
                        help='Validation performance')
    parser.add_argument('-t', '--test-pred',
                        type=str, required=True,
                        help='File to save predictions for the test set')
    parser.add_argument('-b', '--bag-of-words',
                        type=str, required=True,
                        help='bag of words trained model')
    parser.add_argument('-s', '--sequence-specific',
                        type=str, required=True,
                        help='sequence specific trained model')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Predictions for every sequence')
    args = parser.parse_args()
    return vars(args)


def compute_nterm_feature_matrix(sequences, split, dinuc=False):
    """Compute the feature matrix"""
    if 0 < split < 23:
        X = degron_pred.binned_bag_of_words(sequences.str[:split],
                                            int(split), n=int(split),
                                            dinuc=dinuc, cterm=False)
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
                                            dinuc=dinuc, cterm=False)
    return X


def main(opts):
    df = pd.read_table(opts['input'])

    # process dataframe
    feature_df = df.drop_duplicates(subset=['Gene_ID']).copy()

    # setup class label
    y = feature_df['Modal Bin'].apply(lambda x: 1 if x<3.5 else 0)

    #######################
    # setup feature matrix for sequence specific
    #######################
    split_vals = [6, 12, 18, 23]
    dinucs = [True, False] #[True, False]
    performance = [['split', 'dinuc', 'auc']]

    for split in split_vals:
        for dinuc in dinucs:
            # sequence-specific training
            X = compute_nterm_feature_matrix(feature_df['Peptide amino acid sequence'],
                                                split=split, dinuc=dinuc)

            # Train / test split
            X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                                train_size=0.7,
                                                                test_size=0.3,
                                                                random_state=101,
                                                                shuffle=True)
            X_train, X_val, y_train, y_val = train_test_split(X_train, y_train,
                                                            train_size=0.7,
                                                            test_size=0.3,
                                                            random_state=101,
                                                            shuffle=True)

            # iterate through deep learning params
            model = RandomForestClassifier(random_state=101, n_estimators=1000)
            model.fit(X_train, y_train)
            prob_val = model.predict_proba(X_val)[:,1]
            score = metrics.roc_auc_score(y_val, prob_val)
            performance.append([split, dinuc, score])

    # compile performance metrics
    performance_df = pd.DataFrame(performance[1:], columns=performance[0])
    performance_df = performance_df.sort_values('auc', ascending=False)
    # save the performance metrics
    performance_df.to_csv(opts['validation_perf'], sep='\t', index=False)

    # extract the params
    top_params = performance_df.iloc[0]
    best_split = top_params['split']
    best_dinuc = top_params['dinuc']

    # train the best model
    X = compute_nterm_feature_matrix(feature_df['Peptide amino acid sequence'],
                                     split=best_split, dinuc=best_dinuc)
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        train_size=0.7, test_size=0.3,
                                                        random_state=101,
                                                        shuffle=True)
    model = RandomForestClassifier(random_state=101, n_estimators=1000)
    model.fit(X_train, y_train)

    # score the test set
    prob_test = model.predict_proba(X_test)[:,1]
    score = metrics.roc_auc_score(y_test, prob_test)
    print('Test ROC AUC = {0}'.format(score))
    # save the test set predictions
    test_df = pd.DataFrame({'prediction': prob_test, 'y': y_test})
    test_df.to_csv(opts['test_pred'], sep='\t', index=False)

    #############################
    # Save models trained on the full data
    #############################
    # sequence specific model
    clf1 = RandomForestClassifier(random_state=101, n_estimators=1000)
    clf1.fit(X, y)
    ypred_clf1 = clf1.predict_proba(X)[:,1]
    feature_df['sequence position specific'] = ypred_clf1
    # bag of words
    X = compute_nterm_feature_matrix(feature_df['Peptide amino acid sequence'],
                                     split=0, dinuc=False)
    clf2 = RandomForestClassifier(random_state=101, n_estimators=1000)
    clf2.fit(X, y)
    ypred_clf2 = clf2.predict_proba(X)[:,1]
    feature_df['bag of words'] = ypred_clf2
    feature_df['regulatory potential'] = feature_df['sequence position specific'] - feature_df['bag of words']
    # save as pickle files
    with open(opts['sequence_specific'], 'wb') as handle:
        pickle.dump(clf1, handle)
    with open(opts['bag_of_words'], 'wb') as handle:
        pickle.dump(clf2, handle)
    # save scores
    feature_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

