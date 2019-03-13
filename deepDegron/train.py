"""
File: train.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Train c-terminal degron model
"""
import pandas as pd
import argparse
import degron_pred
from sklearn.model_selection import train_test_split
import sklearn.metrics as metrics
from sklearn.linear_model import LogisticRegression
import pickle


def parse_arguments():
    info = 'Train a c-terminal degron model'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='training data')
    parser.add_argument('-b', '--bag-of-words',
                        type=str, required=True,
                        help='bag of words trained model')
    parser.add_argument('-s', '--sequence-specific',
                        type=str, required=True,
                        help='sequence specific trained model')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df = pd.read_table(opts['input'])
    feature_df = df.drop_duplicates(subset=['Gene_ID']).copy()

    # setup class label
    y = feature_df['Modal Bin'].apply(lambda x: 1 if x<2.5 else 0)

    #######################
    # setup feature matrix for sequence specific
    #######################
    X = degron_pred.binned_bag_of_words(feature_df['Peptide amino acid sequence'], 23, dinuc=True)
    # Train / test split
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        train_size=0.7, test_size=0.3,
                                                        random_state=101,
                                                        shuffle=True)
    # do simple logistic regression
    clf1 = LogisticRegression(random_state=101)
    clf1.fit(X_train, y_train)
    # predict on test set
    prob_test = clf1.predict_proba(X_test)[:, 1]
    # make a roc curve
    roc_auc = metrics.roc_auc_score(y_test, prob_test)
    print('Sequence-specific ROC AUC={0}'.format(roc_auc))

    #######################
    # setup feature matrix for sequence specific
    #######################
    X = degron_pred.binned_bag_of_words(feature_df['Peptide amino acid sequence'], 1, dinuc=False)
    # Train / test split
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        train_size=0.7, test_size=0.3,
                                                        random_state=101,
                                                        shuffle=True)
    # do simple logistic regression
    clf2 = LogisticRegression(random_state=101)
    clf2.fit(X_train, y_train)
    # predict on test set
    prob_test = clf2.predict_proba(X_test)[:, 1]
    # make a roc curve
    roc_auc = metrics.roc_auc_score(y_test, prob_test)
    print('Bag of words ROC AUC={0}'.format(roc_auc))

    ####################
    # save models
    ###################
    with open(opts['sequence_specific'], 'wb') as handle:
        pickle.dump(clf1, handle)
    with open(opts['bag_of_words'], 'wb') as handle:
        pickle.dump(clf2, handle)

    import IPython ; IPython.embed()

if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

