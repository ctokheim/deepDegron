import varcode
from deepDegron import sequence_context as sc
import collections
import numpy as np
from deepDegron import variants
from deepDegron import utils
from deepDegron import degron_pred

def degron(variant_list,
           tx, deg_intervals,
           nuc_context=3,
           num_simulations=10000):
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    #trinuc_context = [sc.get_chasm_context(nc) for nc in trinuc_context]
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides

    # return if no substitution
    if not var_sub:
        return None

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)

    # figure out how many affect lysines
    num_deg_mut = 0
    for myvariant in var_sub:
        num_deg_mut += utils.overlap_with_intervals(myvariant, deg_intervals)

    # return p-value of 1 if no degrons affected
    if num_deg_mut == 0:
        return 0, 1

    # simulate new mutations
    random_sub = seq_context.random_pos(trinuc_count, num_simulations)

    # evaluate the simulations
    deg_ct_list = []
    tmp_mut_pos = np.hstack([abs_pos for ctxt, cds_pos, abs_pos in random_sub])
    for sim_pos in tmp_mut_pos:
        # get the variant effect for simulated mutations
        sim_variant_subs = variants.get_mutation_info(sim_pos, tx, dna_change_sub)

        # count the degron mutations
        sim_num_mut = 0
        for sim_variant in sim_variant_subs:
            sim_num_mut += utils.overlap_with_intervals(sim_variant, deg_intervals)
        deg_ct_list.append(sim_num_mut)

    # calculate p-value
    deg_pval = np.sum(1 for x in deg_ct_list if x>=num_deg_mut) / float(num_simulations)

    return num_deg_mut, deg_pval


def degron_with_indel(variant_list,
                      tx, deg_intervals,
                      nuc_context=3,
                      num_simulations=10000):
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    #trinuc_context = [sc.get_chasm_context(nc) for nc in trinuc_context]
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides
    var_sub_no_nmd = utils.filter_nmd_subs(var_sub, tx)

    # filter for indel variants
    var_indel = [x for x in variant_list
                 if x.__class__ in utils.indels
                 and x.variant.is_indel]
    dna_change_indel = [[v.variant.contig, v.variant.start, v.variant.ref, v.variant.alt]
                        for v in var_indel]
    var_indel = utils.nmd(var_indel, tx, drop=True)

    # return if no substitution
    if not var_sub and not var_indel:
        return None

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)
    seq_context_indel = sc.SequenceContext(tx, 0)

    # figure out how many affect lysines
    num_deg_mut = 0
    for myvariant in var_sub_no_nmd + var_indel:
        num_deg_mut += utils.overlap_with_intervals(myvariant, deg_intervals)

    # return p-value of 1 if no degrons affected
    if num_deg_mut == 0:
        return 0, 1

    # simulate new mutations
    if var_sub:
        random_sub = seq_context.random_pos(trinuc_count, num_simulations)
        tmp_mut_pos = np.hstack([abs_pos for ctxt, cds_pos, abs_pos in random_sub])
    if var_indel:
        tmp_input = [('None', len(dna_change_indel))]
        random_indel = seq_context_indel.random_pos(tmp_input, num_simulations)
        tmp_mut_pos_indel = random_indel[0][2]

    # evaluate the simulations
    num_sub, num_indel = len(dna_change_sub), len(dna_change_indel)
    degron_ct, iter_sim = 0, 0
    #for sim_pos in tmp_mut_pos:
    for i in range(num_simulations):
        # get the variant effect for simulated mutations
        if var_sub_no_nmd:
            sim_pos = tmp_mut_pos[i,:]
            sim_variant_subs = variants.get_mutation_info(sim_pos, tx, dna_change_sub)
            num_sim_subs = len(sim_variant_subs)
            sim_variant_subs = utils.filter_nmd_subs(sim_variant_subs, tx)
        else:
            num_sim_subs = 0
            sim_variant_subs = []

        # get infor for indel
        if var_indel:
            sim_pos_indel = tmp_mut_pos_indel[i, :]
            sim_variant_indel = variants.get_mutation_info(sim_pos_indel, tx, dna_change_indel)
            num_sim_indels = len(sim_variant_indel)
            sim_variant_indel = utils.nmd(sim_variant_indel, tx, drop=True)
        else:
            num_sim_indels = 0
            sim_variant_indel = []

        # combine together simulated variants
        sim_combined = sim_variant_subs + sim_variant_indel

        if not sim_combined:
            degron_ct += 1
        else:
            # count the degron mutations
            sim_num_mut = 0
            for sim_var in sim_combined:
                sim_num_mut += utils.overlap_with_intervals(sim_var, deg_intervals)
            if sim_num_mut >= num_deg_mut:
                degron_ct += 1

        # update number of iterations
        iter_sim += 1

        # stop if sufficient number of simulations reached
        if degron_ct>=100:
            break

    # compute p-value
    deg_pval = degron_ct / float(iter_sim)

    return num_deg_mut, deg_pval


def site(variant_list,
         tx, ub_intervals,
         nuc_context=3,
         num_simulations=10000):
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    #trinuc_context = [sc.get_chasm_context(nc) for nc in trinuc_context]
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides

    # return if no substitution
    if not var_sub:
        return None

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)

    # figure out how many overlap sites of interest
    num_ub_mut = 0
    for myvariant in var_sub:
        num_ub_mut += utils.overlap_with_intervals(myvariant, ub_intervals)

    # return p-value of 1 if no degrons affected
    if num_ub_mut == 0:
        return 0, 1

    # simulate new mutations
    random_sub = seq_context.random_pos(trinuc_count, num_simulations)

    # evaluate the simulations
    ub_ct, iter_sim = 0, 0
    tmp_mut_pos = np.hstack([abs_pos for ctxt, cds_pos, abs_pos in random_sub])
    for sim_pos in tmp_mut_pos:
        # get the variant effect for simulated mutations
        sim_variant_subs = variants.get_mutation_info(sim_pos, tx, dna_change_sub)

        # count the degron mutations
        sim_num_mut = 0
        for sim_variant in sim_variant_subs:
            sim_num_mut += utils.overlap_with_intervals(sim_variant, ub_intervals)

        # count if exceeds observed statistic
        if sim_num_mut >= num_ub_mut:
            ub_ct += 1

        # update iteration counter
        iter_sim += 1

        # stop if sufficient number of simulations reached
        if ub_ct>=100:
            break

    # calculate p-value
    ub_pval = ub_ct / float(iter_sim)

    return num_ub_mut, ub_pval


def terminal_degron(variant_list,
                    tx, clf1, clf2,
                    model='cterm',
                    nuc_context=1.5,
                    num_simulations=10000):
    """Simulate the effect of mutations on terminal degrons. Handles both c-terminal and and n-terminal degrons.

    """
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    trinuc_context = [sc.get_chasm_context(nc) for nc in trinuc_context]
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides
    var_sub_no_nmd = utils.filter_nmd_subs(var_sub, tx)

    # filter for indel variants
    var_indel = [x for x in variant_list
                 if x.__class__ in utils.indels
                 and x.variant.is_indel]
    dna_change_indel = [[v.variant.contig, v.variant.start, v.variant.ref, v.variant.alt]
                        for v in var_indel]
    var_indel = utils.nmd(var_indel, tx, drop=True)

    # return if no substitution
    if not var_sub and not var_indel:
        return None

    # figure out the affect on cterminal degrons
    delta_prob = degron_pred.delta_prob(var_sub_no_nmd+var_indel, tx, clf1, clf2, model=model)
    seqs, _ = utils.process_var_seq(var_sub_no_nmd+var_indel, model=model)
    num_impactful_muts = len(seqs)

    # skip if no terminal variants
    if not delta_prob:
        return 0, 1, 0

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)
    seq_context_indel = sc.SequenceContext(tx, 0)

    # simulate new mutations
    if var_sub:
        random_sub = seq_context.random_pos(trinuc_count, num_simulations)
        tmp_mut_pos = np.hstack([abs_pos for ctxt, cds_pos, abs_pos in random_sub])
    # there is no sequence context for indels
    if var_indel:
        #tmp_input = [('None', len(var_indel))]
        tmp_input = [('None', len(dna_change_indel))]
        random_indel = seq_context_indel.random_pos(tmp_input, num_simulations)
        tmp_mut_pos_indel = random_indel[0][2]

    # evaluate the simulations
    num_sub, num_indel = len(dna_change_sub), len(dna_change_indel)
    delta_prob_ct, iter_sim = 0, 0
    for i in range(num_simulations):
        # get info for substitutions
        if var_sub_no_nmd:
            sim_pos = tmp_mut_pos[i, :]
            sim_variant_subs = variants.get_mutation_info(sim_pos, tx, dna_change_sub)
            num_sim_subs = len(sim_variant_subs)
            # filter based on nmd
            sim_variant_subs = utils.filter_nmd_subs(sim_variant_subs, tx)
        else:
            num_sim_subs = 0
            sim_variant_subs = []

        # get info for indel
        if var_indel:
            sim_pos_indel = tmp_mut_pos_indel[i, :]
            sim_variant_indel = variants.get_mutation_info(sim_pos_indel, tx, dna_change_indel)
            num_sim_indels = len(sim_variant_indel)
            # filter based on nmd
            sim_variant_indel = utils.nmd(sim_variant_indel, tx, drop=True)
        else:
            num_sim_indels = 0
            sim_variant_indel = []

        # combine together the simulated variants
        sim_combined = sim_variant_subs + sim_variant_indel

        if not sim_combined:
            delta_prob_ct += 1
        else:
            # get scores from simulations
            sim_delta_prob = degron_pred.delta_prob(sim_combined, tx, clf1, clf2, model=model)

            # count if exceeds observed statistic
            adjust_factor = (num_sub + num_indel) / (num_sim_subs + num_sim_indels)
            #if abs(sim_delta_prob*adjust_factor) >= abs(delta_prob):
            if (sim_delta_prob*adjust_factor) <= delta_prob:
                delta_prob_ct += 1

        # update number of iterations
        iter_sim += 1

        # stop if sufficient number of simulations reached
        if delta_prob_ct>=100:
            break

    # compute p-value
    delta_prob_pval = delta_prob_ct / float(iter_sim)

    return delta_prob, delta_prob_pval, num_impactful_muts


def clustered_truncation(variant_list,
                         tx,
                         nuc_context=3,
                         num_simulations=10000,
                         side='lesser'):
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    #trinuc_context = [sc.get_chasm_context(nc) for nc in trinuc_context]
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides
    var_sub_no_nmd = utils.filter_nmd_subs(var_sub, tx)

    # filter for indel variants
    var_indel = [x for x in variant_list
                 if x.__class__ in utils.indels
                 and x.variant.is_indel]
    dna_change_indel = [[v.variant.contig, v.variant.start, v.variant.ref, v.variant.alt]
                        for v in var_indel]
    var_indel = utils.nmd(var_indel, tx, drop=True)

    # return if no substitution
    if not var_sub and not var_indel:
        return None

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)
    seq_context_indel = sc.SequenceContext(tx, 0)

    # figure out how many are trunc muts
    num_trunc = 0
    for myvariant in var_sub_no_nmd + var_indel:
        is_trunc = myvariant.__class__ in utils.trunc_types
        if is_trunc:
            num_trunc += 1

    # return p-value of 1 if no degrons affected
    if num_trunc == 0:
        return 0, 1

    # simulate new mutations
    if var_sub:
        random_sub = seq_context.random_pos(trinuc_count, num_simulations)
        tmp_mut_pos = np.hstack([abs_pos for ctxt, cds_pos, abs_pos in random_sub])
    if var_indel:
        tmp_input = [('None', len(dna_change_indel))]
        random_indel = seq_context_indel.random_pos(tmp_input, num_simulations)
        tmp_mut_pos_indel = random_indel[0][2]

    # evaluate the simulations
    num_sub, num_indel = len(dna_change_sub), len(dna_change_indel)
    trunc_ct, iter_sim = 0, 0
    #for sim_pos in tmp_mut_pos:
    for i in range(num_simulations):
        # get the variant effect for simulated mutations
        if var_sub_no_nmd:
            sim_pos = tmp_mut_pos[i,:]
            sim_variant_subs = variants.get_mutation_info(sim_pos, tx, dna_change_sub)
            num_sim_subs = len(sim_variant_subs)
            sim_variant_subs = utils.filter_nmd_subs(sim_variant_subs, tx)
        else:
            num_sim_subs = 0
            sim_variant_subs = []

        # get infor for indel
        if var_indel:
            sim_pos_indel = tmp_mut_pos_indel[i, :]
            sim_variant_indel = variants.get_mutation_info(sim_pos_indel, tx, dna_change_indel)
            num_sim_indels = len(sim_variant_indel)
            sim_variant_indel = utils.nmd(sim_variant_indel, tx, drop=True)
        else:
            num_sim_indels = 0
            sim_variant_indel = []

        # combine together simulated variants
        sim_combined = sim_variant_subs + sim_variant_indel

        if not sim_combined:
            trunc_ct += 1
        else:
            # count the degron mutations
            sim_num_mut = 0
            for sim_var in sim_combined:
                is_trunc = sim_var.__class__ in utils.trunc_types
                if is_trunc:
                    sim_num_mut += 1
            if side=="greater" and sim_num_mut >= num_trunc:
                trunc_ct += 1
            elif side=="lesser" and sim_num_mut <= num_trunc:
                trunc_ct += 1

        # update number of iterations
        iter_sim += 1

        # stop if sufficient number of simulations reached
        if trunc_ct>=100:
            break

    # compute p-value
    trunc_pval = trunc_ct / float(iter_sim)

    return num_trunc, trunc_pval

