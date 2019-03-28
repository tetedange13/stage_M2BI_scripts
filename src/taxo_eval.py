#!/usr/bin/env python3

"""
Aims to containing functions common to several scripts
"""


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

    set_euk = set()
    for euk_name in dict_euk:
        lineage_euk = taxfoo.get_lineage_as_dict(dict_euk[euk_name])
        euk_taxo_level = lineage_euk[taxonomic_cutoff]
        assert(euk_taxo_level)
        set_euk.add(euk_taxo_level)

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


def lca_last_try(old_set_taxids):
    """
    Try to produce a LCA from a set of taxid, by eliminating
    """
    # old_set_taxids = list(map(int, str_list_taxids.split('&')))
    new_set_taxids = set()
    list_ranks = []
    for taxid in old_set_taxids:
        taxid_name = taxfoo.get_taxid_name(taxid)
        if ('metagenome' in taxid_name or 'uncultured' in taxid_name or 
            'unidentified' in taxid_name or 'phage' in taxid_name.lower() or
            taxid_name == 'synthetic construct' or 
            'virus' in taxid_name.lower()):
            pass
            # print("PB:", taxid_name, taxid)
        else:
            new_set_taxids.add(taxid)

    if not new_set_taxids:
        return 'only_trashes'
    return taxfoo.find_lca(new_set_taxids)


def handle_second(str_list_taxids):
    """
    Take an alignment that has been found to be secondary and try to determine
    its taxonomy (i.e. 1 unique taxid)    
    """
    set_taxid_target = set(map(int, str_list_taxids.split('s')))
    lca = taxfoo.find_lca(set_taxid_target)

    if lca == 1: # If LCA searching fails, LCA==1
        lca_attempt = lca_last_try(set_taxid_target)
        if lca_attempt == 'only_trashes':
            return (lca_attempt, )
        elif lca_attempt == 1:
            return ('unsolved_lca_pb', )
        else:
            return ('lca_last_try', lca_attempt)

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
        return ('notDeterminable', 'FP')
    else:
        taxo_name = lineage[taxonomic_cutoff]
        if taxo_name in set_levels_prok:
            return (taxo_name, "TP")
        return (taxo_name, "FP")


def main_func(csv_index_val, two_col_from_csv, sets_levels, taxonomic_cutoff):
    """
    """
    lineage_val = two_col_from_csv.loc[csv_index_val, "lineage"]
    type_align = two_col_from_csv.loc[csv_index_val, "type_align"]
    if lineage_val.startswith(';'): # Normal
        taxid_to_eval = lineage_val.strip(';')
    else: # Secondary (with 1 unique taxid or more) 
        res_second_handling = handle_second(lineage_val)
        remark_eval = res_second_handling[0]
        if len(res_second_handling) == 1: # Problem (only trashes or unsolved)
            return (csv_index_val, remark_eval, 'FP', remark_eval)
        else:
            taxid_to_eval = res_second_handling[1]

    taxo_name, classif = in_zymo(taxid_to_eval, sets_levels, 
                                        taxonomic_cutoff)
    # return (csv_index_val, taxid_to_eval, taxo_name, classif)

    if lineage_val.startswith(';'):
        return (csv_index_val, taxo_name, classif, type_align)
    return (csv_index_val, taxo_name, classif, remark_eval)
