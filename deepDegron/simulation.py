import varcode
import sequence_context as sc
import collections
import numpy as np
import variants
import utils
import degron_pred

def lysine_mutations(variant_list,
                     tx,
                     nuc_context=3,
                     num_simulations=10000):
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides

    # return if no substitution
    if not var_sub:
        return None

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)

    # figure out how many affect lysines
    new_lys = len(['K' for x in var_sub
                   if x.__class__ is varcode.effects.Substitution and
                   x.aa_alt=='K'
                  ])
    lost_lys = len(['K' for x in var_sub
                    if x.__class__ is varcode.effects.Substitution and
                    x.aa_ref=='K'
                   ])

    # return p-value of 1 if no lysines affected
    if new_lys == 0 and lost_lys == 0:
        return 0, 1, 0, 1

    # simulate new mutations
    random_sub = seq_context.random_pos(trinuc_count, num_simulations)

    # evaluate the simulations
    new_lys_ct, lost_lys_ct, iter_sim = 0, 0, 0
    #lost_lys_ct_list, new_lys_ct_list = [], []
    tmp_mut_pos = np.hstack(abs_pos for ctxt, cds_pos, abs_pos in random_sub)
    for sim_pos in tmp_mut_pos:
        # get the variant effect for simulated mutations
        sim_variant_subs = variants.get_mutation_info(sim_pos, tx, dna_change_sub)

        # count the lysine mutations
        sim_new_lys = len(['K' for x in sim_variant_subs
                           if x.__class__ is varcode.effects.Substitution and
                           x.aa_alt=='K'
                          ])
        sim_lost_lys = len(['K' for x in sim_variant_subs
                            if x.__class__ is varcode.effects.Substitution and
                            x.aa_ref=='K'
                           ])

        # count if exceeds observed statistic
        if sim_new_lys >= new_lys:
            new_lys_ct += 1
        if sim_lost_lys >= lost_lys:
            lost_lys_ct += 1

        # update number of iterations
        iter_sim += 1

        # stop if sufficient number of simulations reached
        if new_lys_ct>=100 and lost_lys_ct>=100:
            break
        #new_lys_ct_list.append(num_new_lys)
        #lost_lys_ct_list.append(num_lost_lys)

    # calculate p-value
    new_lys_pval = new_lys_ct / float(iter_sim)
    lost_lys_pval = lost_lys_ct / float(iter_sim)
    #new_lys_pval = np.sum(1 for x in new_lys_ct_list if x>=new_lys) / float(num_simulations)
    #lost_lys_pval = np.sum(1 for x in lost_lys_ct_list if x>=lost_lys) / float(num_simulations)

    return new_lys, new_lys_pval, lost_lys, lost_lys_pval


def degron(variant_list,
           tx, deg_intervals,
           nuc_context=3,
           num_simulations=1000):
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
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
    tmp_mut_pos = np.hstack(abs_pos for ctxt, cds_pos, abs_pos in random_sub)
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


def ub_site(variant_list,
            tx, ub_intervals,
            nuc_context=3,
            num_simulations=10000):
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides

    # return if no substitution
    if not var_sub:
        return None

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)

    # figure out how many affect lysines
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
    tmp_mut_pos = np.hstack(abs_pos for ctxt, cds_pos, abs_pos in random_sub)
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


def cterm_degron(variant_list,
                 tx, clf1, clf2,
                 nuc_context=3,
                 num_simulations=1000):
    """Simulate the affect of mutations on c-terminal degrons.

    """
    # interpet variant context
    var_sub, dna_change_sub, trinuc_context = sc.get_substitution_trinuc_context(variant_list, tx)
    trinuc_count = collections.Counter(trinuc_context).items() # count the trinucleotides

    # return if no substitution
    if not var_sub:
        return None

    # create sequence context obj
    seq_context = sc.SequenceContext(tx, nuc_context)

    # figure out how many affect lysines
    delta_prob = degron_pred.delta_prob(var_sub, clf1, clf2)

    # skip if no c-terminal variants
    if not delta_prob:
        return 0, 1

    # simulate new mutations
    random_sub = seq_context.random_pos(trinuc_count, num_simulations)

    # evaluate the simulations
    delta_prob_ct, iter_sim = 0, 0
    tmp_mut_pos = np.hstack(abs_pos for ctxt, cds_pos, abs_pos in random_sub)
    for sim_pos in tmp_mut_pos:
        # get the variant effect for simulated mutations
        sim_variant_subs = variants.get_mutation_info(sim_pos, tx, dna_change_sub)

        # get scores from simulations
        sim_delta_prob = degron_pred.delta_prob(sim_variant_subs, clf1, clf2)

        import IPython ; IPython.embed()
        raise

        # count if exceeds observed statistic
        if sim_delta_prob >= delta_prob:
            delta_prob_ct += 1

        # update number of iterations
        iter_sim += 1

        # stop if sufficient number of simulations reached
        if delta_prob_ct>=100:
            break

    # compute p-value
    delta_prob_pval = delta_prob_ct / float(iter_sim)

    return delta_prob, delta_prob_pval
