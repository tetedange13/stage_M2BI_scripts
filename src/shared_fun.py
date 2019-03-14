#!/usr/bin/env python3

"""
Aims to containing functions common to several scripts
"""


import src.ncbi_taxdump_utils as taxo_utils


global taxfoo
taxfoo = taxo_utils.NCBI_TaxonomyFoo()


def generate_sets_zymo(taxonomic_cutoff):
    """
    """
    # Path to dump files:
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    nodes_path = to_dbs + "nt_db/taxo_18feb19/nodes.dmp"
    names_path = to_dbs + "nt_db/taxo_18feb19/names.dmp"
    taxfoo.load_nodes_dmp(nodes_path)
    taxfoo.load_names_dmp(names_path)

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
        

def in_zymo(current_taxid, taxonomic_cutoff, tupl_sets):
    """
    Given the rank of a BLAST hit and a sublist of taxonomic levels, deducted
    from the cutoff (for the taxonomic level),  returns a string in 
    ('TP', 'FP', 'TN', 'FN').
    The name of the species is also requiered, to deal with 'strain' cases
    """
    # print(taxfoo.get_taxid_name(current_taxid), end='\t')
    current_lineage = taxfoo.get_lineage_as_dict(current_taxid)
    if taxonomic_cutoff not in current_lineage.keys():
        return 'FP' # Taxo cutoff not in keys, so False Positive ?

    current_taxo_level = current_lineage[taxonomic_cutoff]
    assert(current_taxo_level)

    set_levels_prok, set_levels_euk = tupl_sets
    if current_taxo_level in set_levels_prok: # TP
        return 'TP'
    elif current_taxo_level in set_levels_euk: # TN
        return 'TN'
    else: # False Positive
        return 'FP'
