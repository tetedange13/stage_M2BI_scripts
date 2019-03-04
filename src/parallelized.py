#!/usr/bin/env python3


"""
Contains all functions used to parallelise processes
"""


from src.shared_fun import handle_strain, in_zymo


# FROM CLUST_TAX.PY:
def rk_search(tupl_enum_taxids, query_rank_func):
    """
    To perform parallel local taxonomic rank search
    """
    idx, taxid = tupl_enum_taxids
    #print(taxid)
    
    return ( "clust_" + str(idx), query_rank_func(int(taxid)) )


def FP_search(tupl_blast_res_and_idx, my_df, taxo_cutoff, taxo_levels):
    """
    """
    idx_arg, blast_res = tupl_blast_res_and_idx
    df_index = "clust_" + str(idx_arg)
    name_df = my_df.loc[df_index, "topHit"]
    rank_df = my_df.loc[df_index, "rank"]
    
    if in_zymo(name_df, rank_df, taxo_cutoff, taxo_levels) == 'FP':
        return ('FP', df_index, blast_res)
    
    return ["NOT_FP"]
    

def taxid_mapping(chunk, set_accession):
    """
    """
    # `chunk` will be a list of CSV rows all with the same name column
    # replace this with your real computation
    to_return = []
    for line in chunk:
        _, acc_nb, taxid, _ = line.rstrip('\n').split('\t')

        if acc_nb in set_accession:
            to_return.append((acc_nb, taxid))
            
    if to_return: #List not empty
        return to_return
