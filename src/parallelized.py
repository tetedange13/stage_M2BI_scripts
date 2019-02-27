#!/usr/bin/env python3


"""
Contains all functions used to parallelise processes
"""


# FROM CLUST_TAX.PY:
def rk_search(tupl_enum_taxids, query_rank_func):
    """
    To perform parallel local taxonomic rank search
    """
    idx, taxid = tupl_enum_taxids
    #print(taxid)
    
    return ( "clust_" + str(idx), query_rank_func(int(taxid)) )


def FP_search(tupl_blast_res_and_idx, my_df, taxo_cutoff):
    """
    """
    idx_arg, blast_res = tupl_blast_res_and_idx
    df_index = "clust_" + str(idx_arg)
    name_df = my_df.loc[df_index, "topHit"]
    rank_df = my_df.loc[df_index, "rank"]
    
    if in_zymo(name_df, rank_df, taxo_cutoff, taxo_levels) == 'FP':
        return ('FP', df_index, blast_res)
    
    return ["NOT_FP"]
    
