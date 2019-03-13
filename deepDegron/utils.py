import varcode
import csv
import variants
import pvalue
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


def nmd(variant_list, tx, drop=False):
    """This function is meant to identify whether an early stop codon
    will cause NMD."""
    # length of normal protein sequence
    normal_prot_len = len(tx.protein_sequence)

    # figure out the length of the last coding exon interval
    last_coding_intvl = tx.coding_sequence_position_ranges[-1]
    last_coding_len = last_coding_intvl[1] - last_coding_intvl[0] + 1

    # rule of thumb for NMD is that stop codons within 50nt of the
    # last exon-exon junction are not efficiently degraded
    nmd_insensitive_len = last_coding_len/3. + 50/3.

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


def keyboard_exit_wrapper(func):
    def wrap(self, timeout=None):
        return func(self, timeout=timeout if timeout is not None else 1e100)
    return wrap


def read_degron_intervals(path):
    """Read in the sites of degrons."""
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)

        degron_dict = {}
        for line in myreader:
            gene = line[0]
            start = int(line[1])
            end = int(line[2])

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
    """Check if variant is in any degron."""
    is_sub = myvariant.__class__ is varcode.effects.Substitution
    if is_sub:
        is_in_intvl = variants.is_overlap(myvariant.aa_mutation_start_offset+1, intervals)
        if is_in_intvl:
            return 1
    return 0


def process_degron_results(output_list):
    """Process the results from the degron enrichment analysis."""
    mycols = ['gene', 'num_degron_muts', 'pvalue']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['qvalue'] = pvalue.bh_fdr(output_df['pvalue'])
    return output_df


def process_ub_results(output_list):
    """Process the results from the degron enrichment analysis."""
    mycols = ['gene', 'num_ub_site_muts', 'pvalue']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['qvalue'] = pvalue.bh_fdr(output_df['pvalue'])
    return output_df


def process_lysine_results(output_list):
    """Process the results from the lysine mutation analysis."""
    mycols = ['gene', 'new_lys_muts', 'new_lys_pvalue',
              'lost_lys_muts', 'lost_lys_pvalue']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['new_lys_qvalue'] = pvalue.bh_fdr(output_df['new_lys_pvalue'])
    output_df['lost_lys_qvalue'] = pvalue.bh_fdr(output_df['lost_lys_pvalue'])
    return output_df


def process_cterm_degron_results(output_list):
    """Process the results from cterminal degron mutation analysis."""
    mycols = ['gene', 'delta_reg_potential', 'pvalue']
    output_df = pd.DataFrame(output_list, columns=mycols)
    output_df['qvalue'] = pvalue.bh_fdr(output_df['pvalue'])
    return output_df

