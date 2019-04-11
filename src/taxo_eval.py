#!/usr/bin/env python3

"""
Aims to containing functions common to several scripts
"""


import pandas as pd
import src.ncbi_taxdump_utils as taxo_utils


global taxfoo
taxfoo = taxo_utils.NCBI_TaxonomyFoo()
want_taxo = taxo_utils.default_want_taxonomy

# Path to dump files:
to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
nodes_path = to_dbs + "nt_db/taxo_18feb19/nodes.dmp"
names_path = to_dbs + "nt_db/taxo_18feb19/names.dmp"
taxfoo.load_nodes_dmp(nodes_path)
taxfoo.load_names_dmp(names_path)


def generate_sets_zymo(taxonomic_cutoff):
    """
    """
    # Define the species contained within the Zymo mock community:
    dict_prok = {'Listeria monocytogenes':1639, 
                 'Bacillus subtilis':1423, 
                 'Staphylococcus aureus':1280, 
                 'Escherichia coli':562, 
                 'Lactobacillus fermentum':1613, 
                 'Enterococcus faecalis':1351,
                 'Pseudomonas aeruginosa':287, 
                 'Salmonella enterica':28901}
    dict_euk = {'Cryptococcus neoformans':5207, 
                'Saccharomyces cerevisiae':4932}

    set_prok = set()
    for prok_name in dict_prok:
        lineage_prok = taxfoo.get_lineage_as_dict(dict_prok[prok_name])
        prok_taxo_level = lineage_prok[taxonomic_cutoff]
        assert(prok_taxo_level)
        set_prok.add(prok_taxo_level)
    del prok_name

    set_euk = set()
    for euk_name in dict_euk:
        lineage_euk = taxfoo.get_lineage_as_dict(dict_euk[euk_name])
        euk_taxo_level = lineage_euk[taxonomic_cutoff]
        assert(euk_taxo_level)
        set_euk.add(euk_taxo_level)
    del euk_name

    return (set_prok, set_euk)


def calc_taxo_shift(arg_taxid, taxonomic_cutoff):
    """
    Calculate the 'taxonomic shift', i.e. the number of taxonomic levels of
    difference, between a taxonomic cutoff and the rank of a taxid
    """
    idx_cutoff = want_taxo.index(taxonomic_cutoff) + 1
    sublist_taxo = want_taxo[0:idx_cutoff]
    shift = len(taxfoo.get_lineage(arg_taxid)) - len(sublist_taxo)

    print("BONJ", taxfoo.get_taxid_name(arg_taxid), shift)
    return shift
    # if shift <= 0:
    #     pass
    
    # print(len(taxfoo.))


def get_majo(list_of_things, cutoff_majo):
    """
    Given a list of things (could be numbers of list of strings), determine
    wether there is one whose frequency is up to a given cutoff (so that can be
    considered as majority)
    """
    val_counts = pd.Series(list_of_things).value_counts()
    freq_counts = val_counts/len(list_of_things)
    are_up_to_cutoff  = freq_counts > cutoff_majo
    nb_majo = sum(are_up_to_cutoff)
    assert(nb_majo < 2)
    # if len(set(val_counts)) == 1:
    #     print([taxfoo.get_taxid_name(taxid) for taxid in val_counts.index])

    to_return = [val for val in val_counts]
    if nb_majo == 1: # 1 unique majoritary
        return (str(freq_counts.idxmax()), to_return)
    return ('noMajo', to_return) # NO majoritary


def majo_voting(list_taxid_target, taxonomic_cutoff):
    """
    
    """
    dict_taxid2ancester = {}
    for taxid in set(list_taxid_target):
        lineage = taxfoo.get_dict_lineage_as_taxids(taxid)
        if taxonomic_cutoff in lineage.keys():
            dict_taxid2ancester[taxid] = lineage[taxonomic_cutoff]
        else:
            # dict_taxid2ancester[taxid] = 'notInKey'
            dict_taxid2ancester[taxid] = 'notInKey_' + str(taxid)
    del taxid

    list_taxids_ancesters = [dict_taxid2ancester[taxid] 
                             for taxid in list_taxid_target]
    majo_ancester, _ = get_majo(list_taxids_ancesters, 0.5)

    if majo_ancester == 'noMajo': # Still NO majoritary
        return ('no_majo_found', )
    else:
        if 'notInKey' in majo_ancester:
            return ('majo_notInKey', majo_ancester.split('_')[1])
        # assert (majo_ancester != 'notInKey')
        return ('majo_found', majo_ancester) # Majo found only after the 2nd round


def make_lca(list_taxid_target):
    """
    Take an alignment that has been found to be secondary and try to determine
    its taxonomy (i.e. 1 unique taxid)    
    """
    set_taxid_target = set(list_taxid_target)
    lca = taxfoo.find_lca(set_taxid_target)

    if lca == 1: # If LCA searching fails, LCA==1
        return ('lca_reached_root', lca)

    return ('lca', lca)


def in_zymo(taxo_taxid, tupl_sets, taxonomic_cutoff):
    """
    Given the taxid of a read (lowest one in the taxo), determine if the 
    organism belongs to the Zymo mock comm, at a given taxonomic cutoff 
    """
    set_levels_prok, set_levels_euk = tupl_sets
    lineage = taxfoo.get_lineage_as_dict(taxo_taxid)
    taxo_levels = lineage.keys()

    if taxonomic_cutoff not in taxo_levels:
        # return ('notDeterminable', 'FP')
        return (taxfoo.get_taxid_name(int(taxo_taxid)), 'FP', 'notInKeys')
    else:
        taxo_name = lineage[taxonomic_cutoff]
        if taxo_name in set_levels_prok:
            return (taxo_name, "TP", 'true_pos')
        return (taxo_name, "FP", 'misassigned')
