#!/usr/bin/env python3


"""
Contains all functions used to parallelise processes
"""


import src.shared_fun as shared


# FROM STAT_TAX.PY:
def handle_suppl(list_align_obj, dict_conv_seqid2taxid):
    """
    Take an alignment that has been found to be supplementary
    """
    set_taxid_target = set()
    for alignment in list_align_obj:
        taxid_target = dict_conv_seqid2taxid[alignment["ref_name"]]
        assert(taxid_target) # If taxid not 'None'
        set_taxid_target.add(taxid_target)

    # Different target name, but corresponding to the same taxid ?
    if len(set_taxid_target) == 1: 
        return alignment["ref_name"] # Any target name is OK
    return False


def SAM_taxo_classif(align_list, conv_seqid2taxid, taxonomic_cutoff, tupl_sets,
                     cutoff_ratio):
    """
    Parallelized taxonomic classification, from a group (by readID) of mapped 
    reads. Filter out of 
    """
    if len(align_list) > 1: # Chimeric alignment
        return ('suppl', )
        # return ('suppl', handle_suppl(align_list, conv_seqid2taxid))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

    else: # Normal linear case
        align_dict = align_list[0]
        mapq, ratio_len = align_dict["mapq"], align_dict["ratio_len"]
        if ratio_len >= cutoff_ratio:
            # ref_name = alignment.reference_name
            # ref_name = alignment.reference_name.split()[0]
            current_taxid = conv_seqid2taxid[align_dict["ref_name"]]
            return (shared.in_zymo(current_taxid, taxonomic_cutoff, tupl_sets), 
                    mapq)
        else:
            return ('ratio', ratio_len)


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
    
    if shared.in_zymo(name_df, rank_df, taxo_cutoff) == 'FP':
        return ('FP', df_index, blast_res)
    
    return ["NOT_FP"]
    

def taxid_mapping(chunk, set_accession):
    """
    `chunk` will be a list of CSV rows all with the same name column

    set_accession need to contain NCBI accession numbers (not GI, but can
    easily be changed), WITHOUT version number (".one_number.second_number") 
    """
    to_return = []
    for line in chunk:
        acc_nb, _, taxid, _ = line.rstrip('\n').split('\t')
        to_search = acc_nb

        if acc_nb.startswith('NZ_'):
            to_search = acc_nb[3: ]

        if to_search in set_accession:
            to_return.append((acc_nb, taxid))

    if to_return:
        return to_return
