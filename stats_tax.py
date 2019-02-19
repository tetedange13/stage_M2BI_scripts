#!/usr/bin/env python3


"""
"""


import sys, os
import pandas as pd


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


def in_zymo(sp_name, sp_rank, taxo_level_cutoff):
    """
    Given the rank of a BLAST hit and a sublist of taxonomic levels, deducted
    from the cutoff (for the taxonomic level),  returns a string in 
    ('TP', 'FP', 'TN', 'FN').
    The name of the species is also requiered, to deal with 'strain' cases
    """
    taxo_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 
                   'genus', 'species', 'subspecies'] + ["strain"]
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
        
        
# MAIN:
if __name__ == "__main__":
    # CALCULATION OF STATISTICS:
    dict_stats = {'TN':0, 'FN':0, 'TP':0, 'FP':0}
    dict_str = {'TN': "Assigned et un de nos Euk",
                'TP': "Assigned et une de nos bactos !",
                'FN': "Not assigned ! (cutoff plus bas dans l'arbre)",
                'FP': "Assigned, mais appartient pas à la Zymo"}
    
    
    rownames_df = df_hits.index
    #print(df_hits) ; sys.exit()
    for clust in rownames_df:
        if "clust" in clust:
            splitted_remarks = df_hits.loc[clust, "remarks"].split()
            print("\n", clust, " | ", "E-VAL:", splitted_remarks[-1],
            " | ", "SCORE:", splitted_remarks[-2])
                 
            print("NAME:", df_hits.loc[clust, "topHit"], " | ", 
                  "RANK:", df_hits.loc[clust, "rank"])

            res = in_zymo(df_hits.loc[clust, "topHit"], 
                          df_hits.loc[clust, "rank"], 
                          taxonomy_level_cutoff)
            

            if res == "FP":
                # Look for alternative hits:
                print("Not within the Zymo --> Look for alternative hit!")

                alt_index = "alt_" + clust.split('_')[1] + "_1"
                if alt_index in rownames_df:
                    res_alt = in_zymo(df_hits.loc[alt_index, "topHit"], 
                                      df_hits.loc[alt_index, "rank"], 
                                      taxonomy_level_cutoff)
                    if res_alt != 'FP':
                        print("FOUND:", df_hits.loc[alt_index, "topHit"], 
                              " | ", "RANK:", df_hits.loc[alt_index, "rank"]) 
                        remarks_alt = df_hits.loc[alt_index, "remarks"]
                        splitted_remarks_alt = remarks_alt.split()
                        print(alt_index, " | ", "E-VAL:", 
                              splitted_remarks_alt[-1], " | ", "SCORE:", 
                              splitted_remarks_alt[-2])        
                    else:
                        print("The alternative is still a FP")
                        
                    dict_stats[res_alt] += 1
                    print(dict_str[res_alt])
                        
                else:
                    print("NO alternative hit available")
                    dict_stats[res] += 1 
                    
            
            else:
                print(dict_str[res])
                dict_stats[res] += 1    
        #print(res)        

        #print('\n', df_hits.loc[clust, :])
    
    print(dict_stats)


