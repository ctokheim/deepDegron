import numpy as np
import scipy.stats as stats
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
import pickle

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
    chunk_size = round(23/splits)
    xlist = []
    for j in range(splits):
        if j == (splits-1):
            tmp_x = vectorizer.transform(pep_sequence.str[:-j*chunk_size-1])
        else:
            tmp_x = vectorizer.transform(pep_sequence.str[-(j+1)*chunk_size-1:-j*chunk_size-1])
        xlist.append(tmp_x.toarray())
    out = np.concatenate(xlist, axis=1)

    # add c-terminal dinucleotides
    if dinuc:
        tmp = np.vstack(pep_sequence.str[-3:-1].apply(kmer_count).values)
        out = np.concatenate([out, tmp], axis=1)

    return out


def load_classifier(file_path):
    with open(file_path, 'rb') as handle:
        clf = pickle.load(handle)
    return clf
