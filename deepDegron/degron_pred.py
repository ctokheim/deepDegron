import numpy as np
import scipy.stats as stats
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
import pickle
import utils

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
        kmer_dict[tmp_dimer] += 1
    output_list = np.array(list(kmer_dict.values()))
    return output_list


def binned_bag_of_words(pep_sequence, splits, n=23, dinuc=False):
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


def compute_feature_matrix(sequences, split, dinuc=False):
    """Compute the feature matrix"""
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


def load_classifier(file_path):
    with open(file_path, 'rb') as handle:
        clf = pickle.load(handle)
    return clf


def delta_prob(variants, tx, clf1, clf2):
    """Calculate the difference between a position specific
    model and a "bag of words" model"""
    # fetch c-terminal sequence
    cterm_seq = []
    for v in variants:
        if v.mutant_protein_sequence:
            if type(v) in utils.indels:
                cterm_seq.append(v.mutant_protein_sequence[-23:])
            elif type(v) in utils.base_substitutions:
                if v.aa_mutation_start_offset>(len(v.transcript.protein_sequence) - 23):
                    cterm_seq.append(v.mutant_protein_sequence[-23:])

    # return None if no variants
    if not cterm_seq:
        return 0
    # return None if U in protein sequence
    if 'U' in tx.protein_sequence[-23:]:
        return 0

    # construct dataframe
    result_df = pd.DataFrame({'seq': cterm_seq})

    # create feature matrix
    X = compute_feature_matrix(result_df['seq'], 6, dinuc=True)
    X2 = compute_feature_matrix(result_df['seq'], 0, dinuc=False)

    # predict scores
    result_df['prob'] = clf1.predict_proba(X)[:, 0]
    result_df['prob2'] = clf2.predict_proba(X2)[:, 0]
    result_df['delta prob'] = result_df['prob'] - result_df['prob2']

    # adjust for baseline score
    wt_seq = tx.protein_sequence[-23:]
    wt_df = pd.DataFrame({'seq': [wt_seq]})
    # create feature matrix
    X = compute_feature_matrix(wt_df['seq'], 6, dinuc=True)
    X2 = compute_feature_matrix(wt_df['seq'], 0, dinuc=False)
    wt_df['prob'] = clf1.predict_proba(X)[:, 0]
    wt_df['prob2'] = clf2.predict_proba(X2)[:, 0]
    wt_df['delta prob'] = wt_df['prob'] - wt_df['prob2']
    baseline = wt_df['delta prob'].iloc[0]

    # add up scores
    delta_prob_sum = (result_df['delta prob'] - baseline).sum()

    return delta_prob_sum
