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


def taxo_from_taxid(arg_taxid):
    """
    Given a taxid, return a str of the complete taxonomy, with the GreenGenes
    format:
    k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__agalactiae
    """
    lineage = taxfoo.get_lineage_as_dict(arg_taxid)
    possible_lvls = lineage.keys()
    list_lvls = []

    for taxo_lvl in want_taxo:
        if taxo_lvl in possible_lvls:
            first_letter = taxo_lvl[0]
            if taxo_lvl == "superkingdom":
                first_letter = 'k'

            taxo_str = '-'.join(lineage[taxo_lvl].split())

            if taxo_lvl == 'species':
                found_taxa_lvls = set(list_lvls)
                # if len(found_taxa_lvls) == 3 and 'Other' in set(list_lvls):
                if 'Other' in found_taxa_lvls: # If taxo is found incomplete
                    pass # taxo_str kept as a jointure of all 'species' words
                else: # Else we keep only all 'species' words except the 1st one
                    taxo_str = '-'.join(lineage[taxo_lvl].split()[1:])

            list_lvls.append(first_letter +'__' + taxo_str)

        else:
            list_lvls.append('Other')

    # DO NOT add 'Root', cuz pb of length with summarize_taxa.py
    return ';'.join(list_lvls)


def generate_df_zymo():
    """
    """
    # Define the species contained within the Zymo mock community:
    # dict_name2taxid = {'Listeria monocytogenes':1639, 
    #                    'Bacillus subtilis':1423, 
    #                    'Staphylococcus aureus':1280, 
    #                    'Escherichia coli':562, 
    #                    'Lactobacillus fermentum':1613, 
    #                    'Enterococcus faecalis':1351,
    #                    'Pseudomonas aeruginosa':287, 
    #                    'Salmonella enterica':28901}
    # dict_euk = {'Cryptococcus neoformans':5207, 
    #             'Saccharomyces cerevisiae':4932}

    dict_prok = {1639 : 15.9, 
                 1423 : 15.7, 
                 1280 : 13.3, 
                 562 : 10.0, 
                 1613 : 18.8, 
                 1351 : 10.4,
                 287 : 4.6, 
                 28901 : 11.3}
    set_taxids_prok = set(dict_prok.keys())
    dict_tmp = {}

    for taxid in set_taxids_prok:
        dict_tmp[taxid] = taxo_from_taxid(taxid).split(';')
        dict_tmp[taxid] += [dict_prok[taxid]]
        # print(list_taxa_names)
    del taxid

    return pd.DataFrame.from_dict(data=dict_tmp, 
                                  columns=want_taxo+['rel_count'], 
                                  orient='index')


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
    # print(val_counts)
    return ('noMajo', to_return) # NO majoritary


def get_majo_lca(list_of_things, cutoff_discard):
    """
    Given a list of taxids, count occurences of each taxid, discard the taxids
    that have frequecy below a given cutoff and proceed to a LCA seach on the
    remaining taxids
    """
    val_counts = pd.Series(list_of_things).value_counts()
    freq_counts = val_counts/len(list_of_things)
    are_up_to_cutoff  = freq_counts > (1 - cutoff_discard)
    nb_majo = sum(are_up_to_cutoff)
    to_return = [val for val in val_counts]
    if nb_majo > 1:
        print('NB_MAJO:', nb_majo)
    if nb_majo == 0:
        return ('noMajo', to_return) # NO majoritary
    elif nb_majo == 1: # 1 unique majoritary
        return (freq_counts.idxmax(), to_return)
    else:
        return taxfoo.find_lca(freq_counts[are_up_to_cutoff].index, to_return)


def majo_voting(list_taxid_target, taxonomic_cutoff):
    """
    Proceed to the majo determination for a list of given taxids, at the a
    given taxonomic cutoff
    """
    dict_taxid2ancester = {}
    for taxid in set(list_taxid_target):
        lineage = taxfoo.get_dict_lineage_as_taxids(taxid)
        if 'species' in lineage.keys():
            dict_taxid2ancester[taxid] = lineage['species']
        else:
            dict_+taxid2ancester[taxid] = taxid
            if taxfoo.get_taxid_rank(taxid) == 'no rank':
                print(taxfoo.get_taxid_name(taxid))
    del taxid

    list_taxids_ancesters = [dict_taxid2ancester[taxid] 
                             for taxid in list_taxid_target]
    # majo_ancester, _ = get_majo(list_taxids_ancesters, 0.5)
    majo_ancester, _ = get_majo_lca(list_taxids_ancesters, 0.1)

    if majo_ancester == 'noMajo': # Still NO majoritary
        return ('no_majo_found', )
    else:
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


def in_zymo(taxo_taxid, set_levels_prok, taxonomic_cutoff):
    """
    Given the taxid of a read (lowest one in the taxo), determine if the 
    organism belongs to the Zymo mock comm, at a given taxonomic cutoff 
    """
    lineage = taxfoo.get_lineage_as_dict(taxo_taxid)

    if not lineage: # "cannot find taxid {a_taxid}; quitting." --> empty dict
        return ('notDeterminable', 'FP', 'taxid_unknown')
    else:
        taxo_levels = lineage.keys()

        if taxonomic_cutoff not in taxo_levels:
            # return ('notDeterminable', 'FP')
            return (taxfoo.get_taxid_name(int(taxo_taxid)), 'FP', 'notInKeys')
        else:
            taxo_name = lineage[taxonomic_cutoff]
            if taxo_name in set_levels_prok:
                return (taxo_name, "TP", 'true_pos')
            return (taxo_name, "FP", 'misassigned')
