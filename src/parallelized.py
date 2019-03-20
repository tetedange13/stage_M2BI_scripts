#!/usr/bin/env python3


"""
Contains all functions used to parallelise processes
"""


import src.shared_fun as shared


# FROM STAT_TAX.PY:
def handle_second(tupl_dict_item, dict_conv_seqid2taxid):
    """
    Take an alignment that has been found to be secondary
    """
    readID, list_align_obj = tupl_dict_item

    # Check if they are all secondaries (1st one can be the representative):
    assert(all(map(lambda dico: not dico["is_suppl"], list_align_obj)))
    assert(all(map(lambda dico: dico["is_second"], list_align_obj[1: ])))
    # secondaries = [dico["is_second"] for dico in list_align_obj]
    # if list_align_obj[0]["is_second"]:
    #     assert(all(secondaries))
    # else:
    #     assert(all(secondaries[1: ]))

    # if len(list_align_obj) != 6:
    print("SALUT:", len(list_align_obj))
    mapq = list_align_obj[0]["mapq"]
    if mapq > 1:
        print(readID, mapq, [align_obj["AS"] for align_obj in list_align_obj])
        
    list_taxid_target = []
    for alignment in list_align_obj:
        taxid_target = dict_conv_seqid2taxid[alignment["ref_name"]]
        assert(taxid_target) # If taxid not 'None'
        list_taxid_target.append(taxid_target)

    # Different target name, but corresponding to the same taxid ?
    set_taxid_target = set(list_taxid_target)
    if len(set_taxid_target) == 1: # Mergeable directly
        # print("MERGE !")
        return (next(iter(set_taxid_target)), mapq) # Return this unique taxid
    else:
        lca = shared.taxfoo.find_lca(set_taxid_target)
        assert(lca != 1)
        return (lca, mapq)


def SAM_taxo_classif(tupl_dict_item, conv_seqid2taxid, taxonomic_cutoff, 
                     tupl_sets, cutoff_ratio):
    """
    Parallelized taxonomic classification, from a group (by readID) of mapped 
    reads. Filter out of 
    """
    readID, align_list = tupl_dict_item
    nb_alignments_for_readID = len(align_list)
    if nb_alignments_for_readID > 1: # Secondary alignment
        current_taxid, mapq = handle_second(tupl_dict_item, conv_seqid2taxid)
        type_alignment = 'second'
        # return ('second', )
        # return ('second', handle_second(tupl_dict_item, conv_seqid2taxid),
        #         nb_alignments_for_readID)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

    else: # Normal linear case
        align_dict = align_list[0]
        mapq, ratio_len = align_dict["mapq"], align_dict["ratio_len"]
        if ratio_len >= cutoff_ratio:
            current_taxid = conv_seqid2taxid[align_dict["ref_name"]]
            type_alignment = 'linear'
        else:
            return (readID, 'ratio', str(ratio_len), mapq)

    # lineage = shared.taxfoo.get_lineage_as_dict(current_taxid)
    lineage_as_taxids = shared.taxfoo.get_lineage_as_taxids(current_taxid)
    # if taxonomic_cutoff in lineage:
    #     return type_alignment + [lineage[taxonomic_cutoff]] + [mapq]
    # return ['FP'] + [mapq] # ??
    return (readID, type_alignment, lineage_as_taxids, mapq)


# FROM CLUST_TAX.PY:
def rk_search(tupl_enum_taxids, query_rank_func):
    """
    To perform parallel local taxonomic rank search
    """
    idx, taxid = tupl_enum_taxids
    #print(taxid)
    
    return ( "clust_" + str(idx), query_rank_func(int(taxid)) )
    

def taxid_mapping(chunk, set_accession):
    """
    `chunk` will be a list of CSV rows all with the same name column

    set_accession need to contain NCBI accession numbers (not GI, but can
    easily be changed), WITHOUT version number (".one_number.second_number") 
    """
    to_return = []
    for line in chunk:
        acc_nb, _, taxid, gi = line.rstrip('\n').split('\t')
        to_search = acc_nb

        if acc_nb.startswith('NZ_'):
            to_search = acc_nb[3: ]

        if to_search in set_accession:
            to_return.append((to_search, taxid))

    if to_return:
        return to_return
