import varcode
from varcode.effects import *
from pyensembl import ensembl_grch37, ensembl_grch38, EnsemblRelease
import csv
import itertools as it
from deepDegron import utils

def read_variant(var_info, tx, ensemb_ref=ensembl_grch37, catch_errors=False):
    """Annotate a variant on a specific transcript."""
    myvar = varcode.Variant(contig=var_info[0],
                            start=int(var_info[1]),
                            ref=var_info[2],
                            alt=var_info[3],
                            ensembl=ensemb_ref)
    try:
        myeffect = predict_variant_effect_on_transcript(myvar, tx)
    except ValueError:
        myeffect = None
        if not catch_errors:
            raise
    return myeffect


def read_maf(file_path, chrom=None, release=75, include_samples=False):
    """Read a MAF file using varcode."""
    # varcode ensembl release
    data = EnsemblRelease(release)
    if release>75:
        ensemb_ref = ensembl_grch38
    else:
        ensemb_ref = ensembl_grch37

    # read in maf from file
    with open(file_path) as handle:
        myreader = csv.reader(handle, delimiter='\t')

        # parse header
        header = next(myreader)
        chrom_ix = header.index('Chromosome')
        start_ix = header.index('Start_Position')
        ref_ix = header.index('Reference_Allele')
        alt_ix = header.index('Tumor_Seq_Allele2')
        tx_ix = header.index('Transcript_ID')
        gene_ix = header.index('Hugo_Symbol')
        if include_samples:
            samp_ix = header.index('Tumor_Sample_Barcode')

        # read data
        maf_list = list(myreader)

    # filter by chromosome, if specified
    if chrom:
        maf_list = [r for r in maf_list if r[chrom_ix]==chrom]

    # sort the list
    key_func = lambda x: (x[gene_ix], x[tx_ix])
    maf_list = sorted(maf_list, key=key_func)

    # add variants to list
    var_read_errors = 0
    variant_dict = {}
    for (g, t), rows in it.groupby(maf_list, key_func):
        # skip if no transcript
        if t == '':
            continue
        try:
            varcode_tx = data.transcript_by_id(t)
        except ValueError:
            print('{} not found!'.format(t))
            continue

        # set defaults
        variant_dict.setdefault(g, {})
        variant_dict[g].setdefault('transcript_name', t)
        variant_dict[g].setdefault('transcript', varcode_tx)
        if include_samples:
            variant_dict[g].setdefault('samples', [])

        # parse each variant
        variant_list = []
        sample_list = []
        for row in rows:
            # skip if N in ref/alt of variant
            if 'N' in row[ref_ix] or 'N' in row[alt_ix] or len(row[ref_ix])>10:
                continue

            # parse variant
            tmp_info = [row[chrom_ix], row[start_ix], row[ref_ix], row[alt_ix]]
            myeffect = read_variant(tmp_info, varcode_tx, ensemb_ref, catch_errors=True)
            if myeffect:
                variant_list.append(myeffect)
                if include_samples:
                    sample_list.append(row[samp_ix])
            else:
                var_read_errors += 1

        variant_dict[g]['variants'] = variant_list
        if include_samples:
            variant_dict[g]['samples'] = sample_list

    print('Number of problematic variants: {}'.format(var_read_errors))

    return variant_dict


def get_mutation_info(pos, tx, original_dna_change):
    """Create variant objects for simulated mutations."""
    start_offset = tx.first_start_codon_spliced_offset
    sim_variants = []
    for i in range(len(pos)):
        contig, start, ref, alt = original_dna_change[i]
        new_pos = pos[i]

        # figure out the new reference seq
        # for deletions
        is_del = len(ref)>0 and len(alt)==0
        if is_del:
            len_del = len(ref)
            tmp_offset_pos = tx.spliced_offset(int(new_pos)) - start_offset
            if tx.strand == '+':
                new_ref = tx.coding_sequence[tmp_offset_pos:tmp_offset_pos+len_del]
            elif tx.strand == '-':
                new_ref = tx.coding_sequence[tmp_offset_pos:tmp_offset_pos+len_del]
                new_ref = ''.join(utils.base_pairing[x]
                                  for x in new_ref)[::-1]
            new_var_info = [contig, new_pos, new_ref, alt]
        else:
            new_var_info = [contig, new_pos, ref, alt]

        try:
            sim_variants.append(read_variant(new_var_info, tx))
        except:
            pass
    return sim_variants


def is_overlap(var_pos, intervals):
    """Overlap a position with a set of intervals"""
    for s, e in intervals:
        if s<=var_pos<=e: return True
    return False
