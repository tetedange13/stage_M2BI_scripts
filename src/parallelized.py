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

    mapq = list_align_obj[0]["mapq"]
    if False:
        assert(len(set([align_obj["len_align"] for align_obj in list_align_obj])) == 1)
        # print(readID, [align_obj["AS"] for align_obj in list_align_obj],
        #       [align_obj["de"] for align_obj in list_align_obj])
        
    list_taxid_target = []
    for alignment in list_align_obj:
        taxid_target = dict_conv_seqid2taxid[alignment["ref_name"]]
        assert(taxid_target) # If taxid not 'None'
        list_taxid_target.append(taxid_target)

    # Different target name, but corresponding to the same taxid ?
    set_taxid_target = set(list_taxid_target)
    # print("SALUT:", len(set_taxid_target))
    if len(set_taxid_target) == 1: # Mergeable directly
        return ('merged', next(iter(set_taxid_target))) # Return this unique taxid
    else:
        # print(readID, mapq, [align_obj["AS"] for align_obj in list_align_obj], set_taxid_target)
        lca = shared.taxfoo.find_lca(set_taxid_target)
        if lca == 1: # If LCA searching fails, LCA==1
            return ('pb_lca', 
                    '&'.join(str(taxid) for taxid in set_taxid_target))
        return ('lca', lca)


def SAM_taxo_classif(tupl_dict_item, conv_seqid2taxid, tupl_sets, cutoff_ratio):
    """
    Parallelized taxonomic classification, from a group (by readID) of mapped 
    reads. Filter out of 
    """
    readID, align_list = tupl_dict_item
    nb_alignments_for_readID = len(align_list)
    representative = align_list[0]
    mapq, ratio_len = representative["mapq"], representative["ratio_len"]
    len_align = representative["len_align"]

    if nb_alignments_for_readID > 1: # Secondary alignment
        type_alignment, current_taxid = handle_second(tupl_dict_item, 
                                                            conv_seqid2taxid)
        de_list = [dico["de"] for dico in align_list]
        de = max(de_list)
        # de = sum(de_list) / len(de_list)
    else: # Normal linear case
        de = representative["de"]
        if ratio_len >= cutoff_ratio:
            current_taxid = conv_seqid2taxid[representative["ref_name"]]
            type_alignment = 'linear' # i.e. not secondary
        else:
            return (readID, 'ratio_'+str(ratio_len), 'no', mapq, len_align, de)

    if type_alignment == 'pb_lca':
        lineage_as_taxids = current_taxid
    else:
        lineage_as_taxids = shared.taxfoo.get_lineage_as_taxids(current_taxid)

    return (readID, type_alignment, lineage_as_taxids, mapq, len_align, de)


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
