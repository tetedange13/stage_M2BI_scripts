#!/usr/bin/env python3

"""
Aims to containing functions common to several scripts
"""


def handle_strain(sp_name, rank):
    """
    Deals with 'strain' issues (i.e. organism that have been detected beyond
    the species level taxonomicly = 'no rank')
    Cut the name of the organism to keep only the 'Genus species' form and
    attributes the taxonomic level 'strain' to it  
    """
    splitted_sp_name = sp_name.split()
    
    # Very likely a strain if name longer than 2:
    if rank in ('no rank', 'subspecies') and len(splitted_sp_name) > 2: 
        # We keep only the first 2 words and we change the rank value:
        new_sp_name, new_rank_sp = " ".join(splitted_sp_name[0:2]), "strain"
        
        return (new_sp_name, new_rank_sp)
    
    else:
        return (sp_name, rank)
        

def in_zymo(sp_name, sp_rank, taxo_level_cutoff, list_taxo):
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
    sublist_taxo = list_taxo[list_taxo.index(taxo_level_cutoff): ]
    
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
