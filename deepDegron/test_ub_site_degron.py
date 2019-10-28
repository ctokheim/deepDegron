"""
File: test_ub_site_degron.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Test if UB sites have an appropriate lysine
"""
import pandas as pd
import argparse
from variants import *
import utils


def parse_arguments():
    info = 'Test if UB sites have an appropriate lysine'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutation file')
    parser.add_argument('-u', '--ub-sites',
                        type=str, required=True,
                        help='UB site annotation')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    variant_dict = read_maf(opts['input'], '4')
    ub_intvls = utils.read_ub_sites(opts['ub_sites'])

    output_list = [['gene', 'num_matches', 'num_mismatches', 'out_of_range']]
    for gene in variant_dict:
        # variants for a specific gene
        ensembl_tx_name = variant_dict[gene]['transcript_name']
        tx = variant_dict[gene]['transcript']
        variant_list = variant_dict[gene]['variants']

        if ensembl_tx_name in ub_intvls and tx.protein_sequence:
            ub_pos_list = ub_intvls[ensembl_tx_name]

            mismatch_ct = 0 ; match_ct = 0 ; out_of_range_ct = 0
            for i in range(len(ub_pos_list)):
                pos = ub_pos_list[i][0]

                # make sure not out of range
                if pos > len(tx.protein_sequence):
                    out_of_range_ct += 1
                    continue

                aa = tx.protein_sequence[pos-1]
                if aa == 'K':
                    match_ct += 1
                else:
                    mismatch_ct += 1
            output_list.append([gene, match_ct, mismatch_ct, out_of_range_ct])
        if gene == 'KIT':
            import IPython ; IPython.embed()
            raise

    import IPython ; IPython.embed()
    raise


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

