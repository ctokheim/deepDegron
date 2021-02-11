import varcode
import csv
from deepDegron import variants
from deepDegron import pvalue
import pandas as pd

base_pairing = {'A': 'T',
                'T': 'A',
                'a': 't',
                't': 'a',
                'C': 'G',
                'G': 'C',
                'c': 'g',
                'g': 'c',
                '-': '-',  # some people denote indels with '-'
                '': '',    # other people use empty string for indels
                'n': 'n',
                'N': 'N'}

base_substitutions = [varcode.effects.Substitution, varcode.effects.Silent,
                      varcode.effects.AlternateStartCodon, varcode.effects.PrematureStop,
                      varcode.effects.StartLoss, varcode.effects.StopLoss]
indels = [varcode.effects.Deletion, varcode.effects.Insertion, varcode.effects.FrameShift]
nmd_sub_vars = [varcode.effects.PrematureStop]
trunc_types = [
    # varcode.effects.PrematureStop,
    varcode.effects.FrameShift,
    varcode.effects.PrematureStop, varcode.effects.StopLoss]

def rev_comp(seq):
    """Get reverse complement of sequence.
    rev_comp will maintain the case of the sequence.
    Parameters
    ----------
    seq : str
        nucleotide sequence. valid {a, c, t, g, n}
    Returns
    -------
    rev_comp_seq : str
        reverse complement of sequence
    """
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([base_pairing[s] for s in rev_seq])
    return rev_comp_seq


def get_nmd_insensitive_len(tx):
    """Gets the length of the protein that would be insensitive to NMD."""
    # figure out the length of the last coding exon interval
    last_coding_intvl = tx.coding_sequence_position_ranges[-1]
    last_coding_len = last_coding_intvl[1] - last_coding_intvl[0] + 1

    # rule of thumb for NMD is that stop codons within 50nt of the
    # last exon-exon junction are not efficiently degraded
    nmd_insensitive_len = last_coding_len/3. + 50/3.

    return nmd_insensitive_len


def nmd(variant_list, tx, drop=False):
    """This function is meant to identify whether an early stop codon
    will cause NMD."""
    # length of normal protein sequence
    normal_prot_len = len(tx.protein_sequence)
    nmd_insensitive_len = get_nmd_insensitive_len(tx)

    # iterate over each variant
    nmd_sensitive_list = []
    for var in variant_list:
        if not var.mutant_protein_sequence:
            if drop:
                nmd_sensitive_list.append(var)
            else:
                nmd_sensitive_list.append(0.5)
        elif len(var.mutant_protein_sequence) > (normal_prot_len - nmd_insensitive_len):
            if drop:
                nmd_sensitive_list.append(var)
            else:
                nmd_sensitive_list.append(0)
        else:
            if drop:
                pass
            else:
                nmd_sensitive_list.append(1)

    return nmd_sensitive_list


def filter_nmd_subs(var_list, tx):
    """Wrapper function to filter out nonsense mutations which cause NMD from a
    mixture of base substitution mutations. """
    # separate out nonsense mutations
    var_nonsense = [x for x in var_list
                    if x.__class__ in nmd_sub_vars]
    var_other = [x for x in var_list
                 if x.__class__ not in nmd_sub_vars]

    # drop nonsense mutations which cause nmd
    var_nonsense = nmd(var_nonsense, tx, drop=True)

    return var_other+var_nonsense


def keyboard_exit_wrapper(func):
    def wrap(self, timeout=None):
        return func(self, timeout=timeout if timeout is not None else 1e100)
    return wrap


def read_degron_intervals(path, use_tx=False):
    """Read in the sites of degrons."""
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)

        # figure out index
        gene_ix = header.index('gene')
        start_ix = header.index('start')
        end_ix = header.index('end')
        if use_tx:
            gene_ix = header.index('transcript_id')

        degron_dict = {}
        for line in myreader:
            gene = line[gene_ix]
            start = int(line[start_ix])
            end = int(line[end_ix])

            # add interval
            degron_dict.setdefault(gene, [])
            degron_dict[gene].append((start, end))
    return degron_dict


def read_sites(path, num_flank):
    """Read in the sites of degrons."""
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)
        pos_ix = header.index('MOD_RSD')
        enst_ix = header.index('ENST')

        ub_dict = {}
        for line in myreader:
            enst = line[enst_ix]
            pos = int(line[pos_ix])

            # add interval
            ub_dict.setdefault(enst, [])
            ub_dict[enst].append((pos-num_flank, pos+num_flank))
    return ub_dict


def overlap_with_intervals(myvariant, intervals):
    """Check if variant effects any degron."""
    is_sub = myvariant.__class__ is varcode.effects.Substitution
    is_nonsense = myvariant.__class__ is varcode.effects.PrematureStop
    is_fs = myvariant.__class__ is varcode.effects.FrameShift
    is_indel = (myvariant.__class__ is varcode.effects.Deletion) | (myvariant.__class__ is varcode.effects.Insertion)

    if is_sub or is_indel:
        is_in_intvl = variants.is_overlap(myvariant.aa_mutation_start_offset+1, intervals)
        if is_in_intvl:
            return 1
    elif is_nonsense or is_fs:
        var_pos = myvariant.aa_mutation_start_offset+1
        for s, e in intervals:
            if var_pos<e:
                return 1

    return 0


def process_degron_results(output_list):
    """Process the results from the degron enrichment analysis."""
    mycols = ['gene', 'num_degron_muts', 'pvalue', 'total_muts', 'average position']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['qvalue'] = pvalue.bh_fdr(output_df['pvalue'])
    return output_df


def process_trunc_results(output_list):
    """Process the results from the clustered truncation analysis."""
    mycols = ['gene', 'num_truncations', 'pvalue', 'total_muts', 'average position']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['qvalue'] = pvalue.bh_fdr(output_df['pvalue'])
    return output_df


def process_ub_results(output_list):
    """Process the results from the degron enrichment analysis."""
    mycols = ['gene', 'num_ub_site_muts', 'pvalue', 'total_muts', 'average position']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['qvalue'] = pvalue.bh_fdr(output_df['pvalue'])
    return output_df


def process_terminal_degron_results(output_list):
    """Process the results from cterminal degron mutation analysis."""
    mycols = ['gene', 'delta_degron_potential', 'pvalue', 'num_degron_impactful_muts', 'total_muts', 'average position']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['qvalue'] = pvalue.bh_fdr(output_df['pvalue'])
    return output_df


def fetch_seq(seq, model='cterm'):
    """Slices the protein sequence to extract out the n-terminal or c-terminal protein sequence"""
    if model == 'cterm':
        return seq[-23:]
    else:
        return seq[1:24]


def process_var_seq(variants, model='cterm', is_sum=True):
    """This fetches protein sequences for only those that impact the terminal ends of the protein."""
    term_seq = [] ; vars_considered = []
    for v in variants:
        if v.mutant_protein_sequence:
            if model=='cterm' and type(v) in indels+nmd_sub_vars:
                term_seq.append(fetch_seq(v.mutant_protein_sequence, model=model))
                if not is_sum: vars_considered.append(v)
            elif type(v) in base_substitutions:
                if model=='cterm' and v.aa_mutation_start_offset>(len(v.transcript.protein_sequence) - 23):
                    term_seq.append(fetch_seq(v.mutant_protein_sequence, model=model))
                    if not is_sum: vars_considered.append(v)
                elif model=='nterm' and v.aa_mutation_start_offset<=24:
                    term_seq.append(fetch_seq(v.mutant_protein_sequence, model=model))
                    if not is_sum: vars_considered.append(v)
    return term_seq, vars_considered


def fetch_tx_by_gene(gene_name, ensembl_data):
    tx_list = []
    try:
        tx_ids = ensembl_data.transcript_ids_of_gene_name(gene_name)
    except ValueError:
        return tx_list
    for tx_id in tx_ids:
        tx = ensembl_data.transcript_by_id(tx_id)
        if tx.biotype=='protein_coding' and tx.complete:
            tx_list.append(tx)
    return tx_list

