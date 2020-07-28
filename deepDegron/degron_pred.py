import numpy as np
import scipy.stats as stats
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
import pickle
import deepDegron.utils as utils

# load keras
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import Adam
from keras import backend as K

# ignore warnings
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

# don't log deprecation warnings
import logging
logger = logging.getLogger(__name__)  # module logger
logging.getLogger('tensorflow').setLevel(logging.ERROR)
#K.set_session(
    #K.tf.Session(config=K.tf.ConfigProto(intra_op_parallelism_threads=1,
                                         #inter_op_parallelism_threads=2,
                                         #allow_soft_placement=True,
                                         #device_count = {'CPU': 1}))
#)

vocab = ['A', 'R', 'N', 'D', 'C',
         'E', 'Q', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P',
         'S', 'T', 'W', 'Y', 'V']
vectorizer = CountVectorizer(analyzer='char', lowercase=False, vocabulary=vocab)

from collections import OrderedDict
dimers = []
for a1 in vocab:
    for a2 in vocab:
        dimers.append(a1+a2)

def kmer_count(seq):
    kmer_dict = OrderedDict([[aas, 0] for aas in dimers])
    for i in range(len(seq)-1):
        tmp_dimer = seq[i:i+2]
        if 'U' in tmp_dimer: continue  # skip non-DNA letters
        kmer_dict[tmp_dimer] += 1
    output_list = np.array(list(kmer_dict.values()))
    return output_list


def binned_bag_of_words(pep_sequence, splits, n=23, dinuc=False, cterm=True):
    """This function performs bag of words on separate substrings"""
    # split up the bag of words into separate bins
    if splits > 1:
        chunk_size = round(n/splits)
        xlist = []
        for j in range(splits):
            if j == (splits-1):
                tmp_x = vectorizer.transform(pep_sequence.str[:-j*chunk_size])
            else:
                tmp_x = vectorizer.transform(pep_sequence.str[-(j+1)*chunk_size:-j*chunk_size])
            xlist.append(tmp_x.toarray())
        out = np.concatenate(xlist, axis=1)
    else:
        out = vectorizer.transform(pep_sequence)

    # add c-terminal dinucleotides
    if dinuc:
        if cterm:
            tmp = np.vstack(pep_sequence.str[-2:].apply(kmer_count).values)
        else:
            tmp = np.vstack(pep_sequence.str[:2].apply(kmer_count).values)
        out = np.concatenate([out, tmp], axis=1)

    return out


def nterm_binned_bag_of_words(pep_sequence, splits, n=23, dinuc=False):
    """This function performs bag of words on separate substrings"""
    # split up the bag of words into separate bins
    if splits > 1:
        chunk_size = round(n/splits)
        xlist = []
        for j in range(splits):
            if j == (splits-1):
                tmp_x = vectorizer.transform(pep_sequence.str[:-j*chunk_size])
            else:
                tmp_x = vectorizer.transform(pep_sequence.str[-(j+1)*chunk_size:-j*chunk_size])
            xlist.append(tmp_x.toarray())
        out = np.concatenate(xlist, axis=1)
    else:
        out = vectorizer.transform(pep_sequence)

    # add c-terminal dinucleotides
    if dinuc:
        tmp = np.vstack(pep_sequence.str[-2:].apply(kmer_count).values)
        out = np.concatenate([out, tmp], axis=1)

    return out


def compute_nterm_feature_matrix(sequences, split, dinuc=False):
    """Compute the nterm feature matrix"""
    if 0 < split < 23:
        X = binned_bag_of_words(sequences.str[:split],
                                int(split), n=int(split),
                                dinuc=dinuc, cterm=False)
        X2 = binned_bag_of_words(sequences.str[split:],
                                 1, n=23-int(split),
                                 dinuc=False, cterm=False)
        X = np.hstack([X2.toarray(), X])
    elif split == 0:
        X = binned_bag_of_words(sequences,
                                int(split), n=int(split),
                                dinuc=False, cterm=False)
    elif split == 23:
        X = binned_bag_of_words(sequences,
                                int(split), n=int(split),
                                dinuc=dinuc, cterm=False)
    return X


def compute_cterm_feature_matrix(sequences, split, dinuc=False):
    """Compute the cterm feature matrix"""
    if 0 < split < 23:
        X = binned_bag_of_words(sequences.str[-split:],
                                            int(split), n=int(split),
                                            dinuc=dinuc)
        X2 = binned_bag_of_words(sequences.str[:-split],
                                 1, n=23-int(split), dinuc=False)
        X = np.hstack([X2.toarray(), X])
    elif split == 0:
        X = binned_bag_of_words(sequences,
                                int(split), n=int(split),
                                dinuc=False)
    elif split == 23:
        X = binned_bag_of_words(sequences,
                                int(split), n=int(split),
                                dinuc=dinuc)
    return X


def compute_feature_matrix(sequences, split, dinuc=False, model='cterm'):
    """Compute the feature matrix"""
    if model == 'cterm':
        X = compute_cterm_feature_matrix(sequences, split, dinuc=dinuc)
    else:
        X = compute_nterm_feature_matrix(sequences, split, dinuc=dinuc)
    return X


def train_ff_nn(features, y, size=32, dropout=0.5, layers=2, lr=0.001, epochs=40):
    """Train feed-forward neural network model."""
    # compile model
    features_shape = features.shape
    model = Sequential()
    model.add(Dense(units=size, activation='relu', input_dim=features_shape[1]))
    if dropout != 0:
        model.add(Dropout(dropout))
    for i in range(layers-1):
        model.add(Dense(units=size, activation='relu'))
    model.add(Dense(units=1, activation='sigmoid'))
    model.compile(optimizer=Adam(lr=lr),
                  loss='binary_crossentropy',
                  metrics=['accuracy'])

    # fit model
    model.fit(features, y.values,
              epochs=epochs, batch_size=128)
    return model


def load_classifier(file_path):
    with open(file_path, 'rb') as handle:
        clf = pickle.load(handle)
    return clf


def delta_prob(variants, tx, clf1, clf2, model='cterm', is_sum=True):
    """Calculate the difference between a position specific
    model and a "bag of words" model"""
    # fetch c-terminal sequence
    term_seq, vars_considered = utils.process_var_seq(variants, model=model, is_sum=is_sum)

    # return None if no variants
    if not term_seq:
        if is_sum: return 0
        else: return [], []
    # return None if U in protein sequence
    if 'U' in utils.fetch_seq(tx.protein_sequence, model=model):
        if is_sum: return 0
        else: return [], []

    # construct dataframe
    result_df = pd.DataFrame({'seq': term_seq})

    # create feature matrix
    X = compute_feature_matrix(result_df['seq'], 6, dinuc=True, model=model)
    X2 = compute_feature_matrix(result_df['seq'], 0, dinuc=False, model=model)

    # predict scores
    result_df['prob'] = clf1.predict_proba(X)[:, 0]
    result_df['prob2'] = clf2.predict_proba(X2)[:, 0]
    result_df['delta prob'] = result_df['prob'] - result_df['prob2']

    # adjust for baseline score
    wt_seq = utils.fetch_seq(tx.protein_sequence, model=model)
    wt_df = pd.DataFrame({'seq': [wt_seq]})
    # create feature matrix
    X = compute_feature_matrix(wt_df['seq'], 6, dinuc=True, model=model)
    X2 = compute_feature_matrix(wt_df['seq'], 0, dinuc=False, model=model)
    wt_df['prob'] = clf1.predict_proba(X)[:, 0]
    wt_df['prob2'] = clf2.predict_proba(X2)[:, 0]
    wt_df['delta prob'] = wt_df['prob'] - wt_df['prob2']
    baseline = wt_df['delta prob'].iloc[0]

    # add up scores
    tmp = result_df['delta prob'] - baseline
    if is_sum:
        delta_prob_sum = tmp.sum()
        return delta_prob_sum
    else:
        return vars_considered, tmp


def delta_prob_raw(variants, tx, clf1, clf2, model='cterm', is_sum=True):
    """Calculate the difference between a position specific
    model and a "bag of words" model"""
    # fetch c-terminal sequence
    term_seq = [] ; vars_considered = []
    for v in variants:
        if v.mutant_protein_sequence:
            if model=='cterm' and type(v) in utils.indels+utils.nmd_sub_vars:
                term_seq.append(utils.fetch_seq(v.mutant_protein_sequence, model=model))
                if not is_sum: vars_considered.append(v)
            elif type(v) in utils.base_substitutions:
                if model=='cterm' and v.aa_mutation_start_offset>(len(v.transcript.protein_sequence) - 23):
                    term_seq.append(utils.fetch_seq(v.mutant_protein_sequence, model=model))
                    if not is_sum: vars_considered.append(v)
                elif model=='nterm' and v.aa_mutation_start_offset<=24:
                    term_seq.append(utils.fetch_seq(v.mutant_protein_sequence, model=model))
                    if not is_sum: vars_considered.append(v)

    # return None if no variants
    if not term_seq:
        if is_sum: return 0
        else: return [], [], []
    # return None if U in protein sequence
    if 'U' in utils.fetch_seq(tx.protein_sequence, model=model):
        if is_sum: return 0
        else: return [], [], []

    # construct dataframe
    result_df = pd.DataFrame({'seq': term_seq})

    # create feature matrix
    X = compute_feature_matrix(result_df['seq'], 6, dinuc=True, model=model)
    X2 = compute_feature_matrix(result_df['seq'], 0, dinuc=False, model=model)

    # predict scores
    result_df['prob'] = clf1.predict_proba(X)[:, 0]

    # adjust for baseline score
    wt_seq = utils.fetch_seq(tx.protein_sequence, model=model)
    wt_df = pd.DataFrame({'seq': [wt_seq]})
    # create feature matrix
    X = compute_feature_matrix(wt_df['seq'], 6, dinuc=True, model=model)
    wt_df['prob'] = clf1.predict_proba(X)[:, 0]
    baseline = wt_df['prob'].iloc[0]

    # add up scores
    tmp = result_df['prob'] - baseline
    if is_sum:
        prob_sum = tmp.sum()
        return prob_sum
    else:
        return vars_considered, tmp, result_df['prob']


def delta_prob_wildtype(tx, clf1, clf2, model='cterm', seq=None):
    """Calculate the difference between a position specific
    model and a "bag of words" model"""
    if seq is None:
        # fetch c-terminal sequence
        wt_seq = utils.fetch_seq(tx.protein_sequence, model=model)
    else:
        wt_seq = seq

    # return None if U in protein sequence
    if 'U' in wt_seq:
        return None

    # adjust for baseline score
    wt_df = pd.DataFrame({'seq': [wt_seq]})
    # create feature matrix
    X = compute_feature_matrix(wt_df['seq'], 6, dinuc=True, model=model)
    X2 = compute_feature_matrix(wt_df['seq'], 0, dinuc=False, model=model)
    wt_df['prob'] = clf1.predict_proba(X)[:, 0]
    wt_df['prob2'] = clf2.predict_proba(X2)[:, 0]
    wt_df['delta prob'] = wt_df['prob'] - wt_df['prob2']

    return wt_df

