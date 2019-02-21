#!/usr/bin/env python3


"""
"""


import sys, os
import pandas as pd
from src.shared_fun import in_zymo


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
        
        
# MAIN:
if __name__ == "__main__":
    # ARGUMENTS:
    to_report_csv = "cDNA_run1_sampl50k_memb250.csv"
    taxo_cutoff = "genus"
    df_hits = pd.read_csv(to_report_csv, sep='\t', index_col=0)

    # CALCULATION OF STATISTICS:
    cutoff_e_val = 0.00001 # 10^-5
    dict_stats = {'TN':0, 'FN':0, 'TP':0, 'FP':0}
    dict_str = {'TN': "Assigned et un de nos Euk",
                'TP': "Assigned et une de nos bactos !",
                'FN': "Not assigned ! (cutoff plus bas dans l'arbre)",
                'FP': "Assigned, mais appartient pas à la Zymo"}
    
    rownames_df = df_hits.index
    rows_clust = [rowname for rowname in rownames_df if "clust" in rowname]

    for clust in rows_clust:
        splitted_remarks = df_hits.loc[clust, "remarks"].split()
        e_val  = splitted_remarks[-1]
        score = splitted_remarks[-2]
        
        print("\n", clust, " | ", "E-VAL:", e_val, " | ", "SCORE:", score)   
        print("NAME:", df_hits.loc[clust, "topHit"], " | ", 
              "RANK:", df_hits.loc[clust, "rank"])
        
        if float(e_val) < cutoff_e_val:
            res = in_zymo(df_hits.loc[clust, "topHit"], 
                          df_hits.loc[clust, "rank"], 
                          taxo_cutoff)
            

            if res == "FP":
                # Look for alternative hits:
                print("Not within the Zymo --> Look for alternative hit!")

                alt_index = "alt_" + clust.split('_')[1] + "_1"
                if alt_index in rownames_df:
                    res_alt = in_zymo(df_hits.loc[alt_index, "topHit"], 
                                      df_hits.loc[alt_index, "rank"], 
                                      taxo_cutoff)
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
        
        
        else: # Not assigned because of the e_val cutoff
            print("Not assigned because of the e_val cutoff")
            dict_stats['FN'] += 1     
    
    print(dict_stats)


