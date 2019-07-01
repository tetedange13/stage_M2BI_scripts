#!/usr/bin/env python3

"""
Aims to containing functions common to several scripts
"""


import sys
import pandas as pd
import src.ncbi_taxdump_utils as taxo_utils
import os.path as osp


global taxfoo
taxfoo = taxo_utils.NCBI_TaxonomyFoo()
want_taxo = taxo_utils.default_want_taxonomy

# Path to dump files:
to_dumps = "dump_files/"
if not osp.isdir(to_dumps):
    print("ERROR: Need to create a 'dump_files/' folder and put dump files" + 
          " within it")
    print()
    sys.exit()

nodes_path = osp.join(to_dumps, "nodes.dmp")
names_path = osp.join(to_dumps, "names.dmp")
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
    del taxid

    return pd.DataFrame.from_dict(data=dict_tmp, 
                                  columns=want_taxo+['rel_count'], 
                                  orient='index')


def get_majo_old(list_of_things, cutoff_majo):
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

    to_return = [val for val in val_counts]
    if nb_majo == 1: # 1 unique majoritary
        return (str(freq_counts.idxmax()), to_return)

    return ('noMajo', to_return) # NO majoritary


def majo_voting_old(list_taxid_target, taxonomic_cutoff):
    """
    Proceed to the majo determination for a list of given taxids, at the a
    given taxonomic cutoff
    """   
    dict_taxid2ancester = {}
    for taxid in set(list_taxid_target):
        lineage = taxfoo.get_dict_lineage_as_taxids(taxid)
        if taxonomic_cutoff in lineage.keys():
            dict_taxid2ancester[taxid] = lineage[taxonomic_cutoff]
        else:
            dict_taxid2ancester[taxid] = 'notInKey_' + str(taxid)
    del taxid

    list_taxids_ancesters = [dict_taxid2ancester[taxid] 
                             for taxid in list_taxid_target]
    if len(set(list_taxids_ancesters)) > 1:
        majo_ancester, _ = get_majo_old(list_taxids_ancesters, 0.5)
        if majo_ancester == 'noMajo': # Still NO majoritary
            return ('no_majo_found', )
        else:
            if 'notInKey' in majo_ancester:
                return ('majo_notInKey', majo_ancester.split('_')[1])
            return ('majo_found', majo_ancester) # Majo found only after the 2nd round

    else:
        return ('single_taxid', list_taxid_target[0])


def remove_minor_lca(list_of_things, cutoff_discard):
    """
    Given a list of taxids, count occurences of each taxid, discard the taxids
    that have frequecy below a given cutoff and proceed to a LCA seach on the
    remaining taxids
    """
    val_counts = pd.Series(list_of_things).value_counts()
    freq_counts = val_counts/len(list_of_things)
    if len(freq_counts) > 2 and len(set(freq_counts)) > 1:
        freq_counts.name='freq'
        test = pd.DataFrame(freq_counts.apply(round, args=(2, ))).reset_index()
        test = test.assign(sp_name=test['index'].apply(taxfoo.get_taxid_name),
                           diff_freq=lambda x: round(x.freq - 1/len(freq_counts), 2))
    are_down_to_cutoff  = freq_counts < cutoff_discard
    nb_majo = len(are_down_to_cutoff) - sum(are_down_to_cutoff)
    to_return = [val for val in val_counts]

    if nb_majo == 0:
        print('READ noMajo !')
        return ('noMajo', None) # NO majoritary
    
    elif nb_majo == 1: # 1 unique majoritary
        return ('uniq_majo', freq_counts.idxmax())
    
    else:
        # '-' before a pdSeries of bool invert the Series:
        return ('minors_rm_lca', 
                taxfoo.find_lca(list(freq_counts[-are_down_to_cutoff].index)))


def majo_voting(list_taxid_target):
    """
    Proceed to the majo determination for a list of given taxids
    This version is independant from any taxo cutoff
    """
    dict_taxid2ancester = {}
    for taxid in set(list_taxid_target):
        lineage = taxfoo.get_dict_lineage_as_taxids(taxid)
        if 'species' in lineage.keys():
            dict_taxid2ancester[taxid] = lineage['species']
        else:
            dict_taxid2ancester[taxid] = taxid
    del taxid

    list_taxids_ancesters = [dict_taxid2ancester[taxid] 
                             for taxid in list_taxid_target]
    if len(set(list_taxids_ancesters)) > 1:
        remark_majo, majo_ancester = remove_minor_lca(list_taxids_ancesters, 
                                                      0.2)
        if remark_majo == 'noMajo': # Still NO majoritary
            return ('no_majo_found', )
        else:
            return (remark_majo, majo_ancester) # Majo found only after the 2nd round
    
    else:
        return ('single_taxid', list_taxid_target[0])


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
    lineage = taxfoo.get_dict_lineage_as_taxids(taxo_taxid)

    if not lineage: # "cannot find taxid {a_taxid}; quitting." --> empty dict
        return (pd.np.nan, 'miss', 'taxid_unknown')
    else:
        taxo_levels = lineage.keys()

        if taxonomic_cutoff not in taxo_levels:
            # return ('notDeterminable', 'miss')
            return (pd.np.nan, 'miss', 'notInKeys')
        else:
            taxid_ancester = lineage[taxonomic_cutoff]
            assert(taxid_ancester)
            taxo_name = taxfoo.get_taxid_name(taxid_ancester)
            if taxo_name in set_levels_prok:
                return (taxid_ancester, 'well', 'true_pos')
            return (taxid_ancester, 'miss', 'misassigned')


def in_zymo_new(taxo_taxid, set_levels_prok, taxonomic_cutoff):
    """
    Given the taxid of a read (lowest one in the taxo), determine if the 
    organism belongs to the Zymo mock comm, at a given taxonomic cutoff 
    """
    if taxo_taxid == 'no_majo_found':
        returned_list = ['no_majo_found', 'miss', 'no_majo_found']

    else:
        lineage = taxfoo.get_dict_lineage_as_taxids(taxo_taxid)

        if not lineage: # "cannot find taxid {a_taxid}; quitting." --> empty dict
            returned_list = [pd.np.nan, 'miss', 'taxid_unknown']
        else:
            taxo_levels = lineage.keys()

            if taxonomic_cutoff not in taxo_levels:
                returned_list = [pd.np.nan, 'miss', 'notInKeys']
            else:
                taxid_ancester = lineage[taxonomic_cutoff]
                assert(taxid_ancester)
                taxo_name = taxfoo.get_taxid_name(taxid_ancester)
                if taxo_name in set_levels_prok:
                    returned_list = [taxid_ancester, 'well', 'true_pos']
                else:
                    returned_list = [taxid_ancester, 'miss', 'misassigned']

    return pd.Series(returned_list)
