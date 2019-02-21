#!/usr/bin/env python3

"""
Aims to containing functions common to several scripts
"""

def in_zymo(sp_name, sp_rank, taxo_level_cutoff):
    """
    Given the rank of a BLAST hit and a sublist of taxonomic levels, deducted
    from the cutoff (for the taxonomic level),  returns a string in 
    ('TP', 'FP', 'TN', 'FN').
    The name of the species is also requiered, to deal with 'strain' cases
    """
    # Define the species contained within the Zymo mock community:
    list_prok = ['Listeria monocytogenes', 'Bacillus subtilis', 
                 'Staphylococcus aureus', 'Escherichia coli', 
                 'Lactobacillus fermentum', 'Enterococcus faecalis',
                 'Pseudomonas aeruginosa', 'Salmonella enterica']
    list_euk = ['Cryptococcus neoformans', 'Saccharomyces cerevisiae']
    sublist_taxo = taxo_levels[taxo_levels.index(taxo_level_cutoff): ]
    
    # We handle the issues associated with 'strains'
    new_name, rew_rank = handle_strain(sp_name, sp_rank)
    
    if rew_rank in sublist_taxo: # Assigned
        if new_name in list_prok:
            return 'TP' # True Positive
        elif new_name in list_euk:
            return 'TN' # True Negative
        else:
            return 'FP' # False Positive
    
    else:
        return 'FN'
