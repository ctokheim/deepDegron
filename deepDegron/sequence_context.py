import numpy as np
import deepDegron.utils as utils

def cds_pos_to_genome(tx):
    """Creates a dictionary that maps CDS position to genome location.

    Parameters
    ----------
    tx : vacode.Transcript

    Returns
    -------
    seqpos2genome : dictionary
    """
    # initialize variables
    exons = tx.coding_sequence_position_ranges
    cds_len = len(tx.coding_sequence)
    strand = tx.strand
    seqpos2genome = {}

    # record genome positions for each sequence position
    seq_pos = 0
    for start, end in exons:
        myrange = range(start, end+1)
        if strand=="-": myrange = reversed(myrange)
        for genome_pos in myrange:
            #if strand == '+':
                #seqpos2genome[seq_pos] = genome_pos
            #elif strand == '-':
                #tmp = cds_len - seq_pos - 1
                #seqpos2genome[tmp] = genome_pos
            seqpos2genome[seq_pos] = genome_pos
            seq_pos += 1

    # add the stop codon to the end
    stop1, stop2, stop3 = sorted(tx.stop_codon_positions)
    if strand=='-': stop1, stop2, stop3 = stop3, stop2, stop1
    # catch weird case where stop codon spans an exon
    num_ovlp = len(set(seqpos2genome.values()) & set([stop1, stop2, stop3]))
    seq_pos -= num_ovlp
    # add to dictionary
    seqpos2genome[seq_pos] = stop1
    seqpos2genome[seq_pos+1] = stop2
    seqpos2genome[seq_pos+2] = stop3

    return seqpos2genome


def get_trinuc_context(cds, offsets):
    """Get the tri-nucleotide context for a particular position"""
    trinuc_list = []
    for offset in offsets:
        trinuc = cds[offset-1:offset+2]
        if offset == 0:
            trinuc = cds[0] + cds[offset:offset+2]
        elif offset == (len(cds)-1):
            trinuc = cds[offset-1:offset+1] + cds[offset]
        trinuc_list.append(trinuc)
    return trinuc_list


def get_substitution_trinuc_context(var_list, tx):
    """Get tri-nucleotide variant information for single base substitutions."""
    start_offset = tx.first_start_codon_spliced_offset

    # prepare variant info
    var_sub = [x for x in var_list
               if x.variant.is_snv and x.__class__ in utils.base_substitutions
                  and x.transcript.transcript_id == tx.transcript_id
               ]
    dna_change = [[v.variant.contig, v.variant.start, v.variant.ref, v.variant.alt]
                  for v in var_sub]

    # use the tx.spliced_offset function to figure out where the cds variants are
    cds_offsets = [tx.spliced_offset(x[1]) - start_offset
                   for x in dna_change]
    trinuc_context = get_trinuc_context(tx.coding_sequence, cds_offsets)

    # return the results in sorted order by nucleotide context
    myorder = sorted(range(len(trinuc_context)), key=trinuc_context.__getitem__)
    var_sub = [var_sub[i] for i in myorder]
    dna_change = [dna_change[i] for i in myorder]
    trinuc_context = [trinuc_context[i] for i in myorder]

    return var_sub, dna_change, trinuc_context


def get_all_context_names(context_num):
    """Based on the nucleotide base context number, return
    a list of strings representing each context.

    Parameters
    ----------
    context_num : int
        number representing the amount of nucleotide base context to use.
    Returns
    -------
        a list of strings containing the names of the base contexts
    """
    if context_num == 0:
        return ['None']
    elif context_num == 1:
        return ['A', 'C', 'T', 'G']
    elif context_num == 1.5:
        return ['C*pG', 'CpG*', 'TpC*', 'G*pA',
                'A', 'C', 'T', 'G']
    elif context_num == 2:
        dinucs = list(set(
            [d1+d2
             for d1 in 'ACTG'
             for d2 in 'ACTG']
        ))
        return dinucs
    elif context_num == 3:
        trinucs = list(set(
            [t1+t2+t3
             for t1 in 'ACTG'
             for t2 in 'ACTG'
             for t3 in 'ACTG']
        ))
        return trinucs


def get_chasm_context(tri_nuc):
    """Returns the mutation context acording to CHASM.
    For more information about CHASM's mutation context, look
    at http://wiki.chasmsoftware.org/index.php/CHASM_Overview.
    Essentially CHASM uses a few specified di-nucleotide contexts
    followed by single nucleotide context.

    Parameters
    ----------
    tri_nuc : str
        three nucleotide string with mutated base in the middle.
    Returns
    -------
    chasm context : str
        a string representing the context used in CHASM
    """
    # check if string is correct length
    if len(tri_nuc) != 3:
        raise ValueError('Chasm context requires a three nucleotide string '
                         '(Provided: "{0}")'.format(tri_nuc))

    # try dinuc context if found
    if tri_nuc[1:] == 'CG':
        return 'C*pG'
    elif tri_nuc[:2] == 'CG':
        return 'CpG*'
    elif tri_nuc[:2] == 'TC':
        return 'TpC*'
    elif tri_nuc[1:] == 'GA':
        return 'G*pA'
    else:
        # just return single nuc context
        return tri_nuc[1]


class SequenceContext(object):
    """The SequenceContext class allows for deciphering sequence context
    and for randomly permuting mutation positions while respecting sequence context.
    """

    def __init__(self, gene_seq, nuc_context, seed=101):
        # handle nuc contexts
        self._init_context(gene_seq, nuc_context)
        self.seed = seed  # seed for random number generator
        context_names = get_all_context_names(nuc_context)
        self.prng_dict = {
            c: np.random.RandomState(seed=self.seed)
            for c in context_names
        }

        # init mapping of cds position to genome
        self.cdspos2genome = cds_pos_to_genome(gene_seq)
        # convert coordinates from CDS positions to genome on a numpy array
        convert_coord = lambda x: self.cdspos2genome[x]
        self.vec_cdspos2genome = np.vectorize(convert_coord)

    def _init_context(self, gene_seq, nuc_context):
        """Initializes attributes defining mutation contexts and their position.

        The self.context2pos and self.pos2context dictionaries map from
        sequence context to sequence position and sequence position to
        sequence context, respectively. These attributes allow for randomly
        sampling of mutation positions while respecting sequence context in the
        randomization-based test.

        Parameters
        ----------
        gene_seq : GeneSequence
            GeneSequence object from the gene_sequence module
        """
        self.context2pos, self.pos2context = {}, {}
        gene_len = len(gene_seq.coding_sequence)  # get length of CDS
        #five_ss_len = 2*len(gene_seq.five_prime_seq)  # total length of 5' splice sites
        #three_ss_len = 2*len(gene_seq.three_prime_seq)  # total length of 3' splice sites

        if nuc_context in [1, 2]:
            # case where context matters
            index_context = int(nuc_context) - 1  # subtract 1 since python is zero-based index
            for i in range(index_context, gene_len):
                nucs = gene_seq.coding_sequence[i-index_context:i+1]
                self.context2pos.setdefault(nucs, [])
                self.context2pos[nucs].append(i)
                self.pos2context[i] = nucs

            # sequence context for five prime splice site
            #for i, five_ss in enumerate(gene_seq.five_prime_seq):
                #first_nucs = five_ss[1-index_context:1+1]
                #second_nucs = five_ss[2-index_context:2+1]
                #first_pos = 2*i + gene_len
                #second_pos = 2*i + gene_len + 1
                #self.context2pos.setdefault(first_nucs, [])
                #self.context2pos[first_nucs].append(first_pos)
                #self.context2pos.setdefault(second_nucs, [])
                #self.context2pos[second_nucs].append(second_pos)
                #self.pos2context[first_pos] = first_nucs
                #self.pos2context[second_pos] = second_nucs
            # sequence context for three prime splice site
            #for i, three_ss in enumerate(gene_seq.three_prime_seq):
                #first_nucs = three_ss[1-index_context:1+1]
                #second_nucs = three_ss[2-index_context:2+1]
                #first_pos = 2*i + gene_len + five_ss_len
                #second_pos = 2*i + gene_len + five_ss_len + 1
                #self.context2pos.setdefault(first_nucs, [])
                #self.context2pos[first_nucs].append(first_pos)
                #self.context2pos.setdefault(second_nucs, [])
                #self.context2pos[second_nucs].append(second_pos)
                #self.pos2context[first_pos] = first_nucs
                #self.pos2context[second_pos] = second_nucs

            # hack solution for context for first nuc
            if gene_seq.coding_sequence and nuc_context > 1:
                self.pos2context[0] = gene_seq.coding_sequence[0] * 2
                self.context2pos.setdefault(gene_seq.coding_sequence[0]*2, [])
                self.context2pos[gene_seq.coding_sequence[0]*2].append(0)
        elif nuc_context in [1.5, 3]:
            # use the nucleotide context from chasm if nuc
            # context is 1.5 otherwise always use a three
            # nucleotide context
            ncontext = nuc_context
            for i in range(1, len(gene_seq.coding_sequence)-1):
                nucs = gene_seq.coding_sequence[i-1:i+2]
                if ncontext == 1.5:
                    context = get_chasm_context(nucs)
                else:
                    context = nucs
                self.context2pos.setdefault(context, [])
                self.context2pos[context].append(i)
                self.pos2context[i] = context

            # sequence context for five prime splice site
            #for i, five_ss in enumerate(gene_seq.five_prime_seq):
                #first_nucs = five_ss[:3]
                #second_nucs = five_ss[1:4]
                #first_pos = 2*i + gene_len
                #second_pos = 2*i + gene_len + 1
                #if ncontext == 1.5:
                    #first_context = get_chasm_context(first_nucs)
                    #second_context = get_chasm_context(second_nucs)
                #else:
                    #first_context = first_nucs
                    #second_context = second_nucs
                #self.context2pos.setdefault(first_context, [])
                #self.context2pos[first_context].append(first_pos)
                #self.context2pos.setdefault(second_context, [])
                #self.context2pos[second_context].append(second_pos)
                #self.pos2context[first_pos] = first_context
                #self.pos2context[second_pos] = second_context
            # sequence context for three prime splice site
            #for i, three_ss in enumerate(gene_seq.three_prime_seq):
                #first_nucs = three_ss[:3]
                #second_nucs = three_ss[1:4]
                #first_pos = 2*i + gene_len + five_ss_len
                #second_pos = 2*i + gene_len + five_ss_len + 1
                #if ncontext == 1.5:
                    #first_context = get_chasm_context(first_nucs)
                    #second_context = get_chasm_context(second_nucs)
                #else:
                    #first_context = first_nucs
                    #second_context = second_nucs
                #self.context2pos.setdefault(first_context, [])
                #self.context2pos[first_context].append(first_pos)
                #self.context2pos.setdefault(second_context, [])
                #self.context2pos[second_context].append(second_pos)
                #self.pos2context[first_pos] = first_context
                #self.pos2context[second_pos] = second_context

            # hack solution for context for first nuc
            if gene_seq.coding_sequence:
                first_nuc = gene_seq.coding_sequence[0] + gene_seq.coding_sequence[:2]
                if ncontext == 1.5:
                    first_context = get_chasm_context(first_nuc)
                else:
                    first_context = first_nuc
                self.pos2context[0] = first_context
                self.context2pos.setdefault(first_context, [])
                self.context2pos[first_context].append(0)
                last_nuc = gene_seq.coding_sequence[-2:] + gene_seq.coding_sequence[-1]
                if ncontext == 1.5:
                    last_context = get_chasm_context(last_nuc)
                else:
                    last_context = last_nuc
                last_pos = len(gene_seq.coding_sequence) - 1
                self.pos2context[last_pos] = first_context
                self.context2pos.setdefault(last_context, [])
                self.context2pos[last_context].append(last_pos)
        else:
            # case where there is no context,
            # mutations occur with uniform probability at each
            # position
            #for i in range(gene_len + five_ss_len + three_ss_len):
            for i in range(gene_len):
                self.pos2context[i] = 'None'
            self.context2pos['None'] = range(gene_len)

    def is_valid_context(self, ctxt):
        """Checks if provided context is valid (previously seen).

        Parameters
        ----------
        ctxt : str
            mutation context
        """
        return ctxt in self.context2pos

    def random_context_pos(self, num, num_permutations, context):
        """Samples with replacement available positions matching the
        sequence context.

        Note: this method does random sampling only for an individual
        sequence context.

        Parameters
        ----------
        num : int
            Number of positions to sample for each permutation. This
            is the number of actually observed mutations having the
            matching sequence context for this gene.
        num_permutations : int
            Number of permutations for permutation test.
        context : str
            Sequence context.

        Returns
        -------
        random_pos : np.array
            num_permutations X num sized array that represents the
            randomly sampled positions for a specific context.
        """
        # make sure provide context is valid
        if not self.is_valid_context(context):
            error_msg = 'Context ({0}) was never seen in sequence.'.format(context)
            raise ValueError(error_msg)

        # make sure sampling is a positive integer
        if num < 1:
            error_msg = ('There must be at least one sample (specified {0}) '
                         'for a context'.format(num))
            raise ValueError(error_msg)

        # randomly select from available positions that fit the specified context
        available_pos = self.context2pos[context]
        random_pos = self.prng_dict[context].choice(available_pos, (num_permutations, num))
        return random_pos

    def random_pos(self, context_iterable, num_permutations):
        """Obtains random positions w/ replacement which match sequence context.

        Parameters
        ----------
        context_iterable: iterable containing two element tuple
            Records number of mutations in each context. context_iterable
            should be something like [('AA', 5), ...].
        num_permutations : int
            Number of permutations used in the permutation test.

        Returns
        -------
        position_list : list
            Contains context string and the randomly chosen positions
            for that context.
        """
        position_list = []
        for contxt, n in context_iterable:
            pos_array = self.random_context_pos(n, num_permutations, contxt)
            genome_pos_array = self.vec_cdspos2genome(pos_array)
            position_list.append([contxt, pos_array, genome_pos_array])
        return position_list
