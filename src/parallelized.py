#!/usr/bin/env python3


"""
Contains all functions used to parallelise processes
"""


import src.shared_fun as shared


# FROM STAT_TAX.PY:
def is_trash(taxid_to_eval):
    """
    """
    taxid_name = shared.taxfoo.get_taxid_name(taxid_to_eval)
    if ('metagenome' in taxid_name or 'uncultured' in taxid_name or 
        'unidentified' in taxid_name or 'phage' in taxid_name.lower() or
        taxid_name == 'synthetic construct' or 
        'virus' in taxid_name.lower()):
        return True
    return False


def SAM_to_CSV(tupl_dict_item, conv_seqid2taxid):
    """
    Prepare the writting of a CSV file summarizing the SAM
    """
    readID, align_list = tupl_dict_item
    nb_alignments_for_readID = len(align_list)
    representative = align_list[0]
    mapq, ratio_len = representative["mapq"], representative["ratio_len"]
    len_align = representative["len_align"]
    common_to_return = [mapq, len_align, ratio_len]

    if nb_alignments_for_readID > 1: # Secondary alignment
        assert(not representative["is_suppl"])
        assert(not representative["is_second"])

        no_suppl_list = [align_obj for align_obj in align_list 
                         if not align_obj["has_SA"]]
        # Check if they are all secondaries (1st one can be the representative):
        # assert(all(map(lambda dico: not dico["is_suppl"], no_suppl_list)))
        # assert(all(map(lambda dico: dico["is_second"], no_suppl_list[1: ])))

        if not no_suppl_list: # Only supplementaries entries for this read
            return (readID, 'only_suppl', 'no')
        # Need to propagate max value of MAPQ to all secondaries:
        if representative["has_SA"]:
            # print(readID, "suppl_as_repr")
            max_mapq = max([align_obj["mapq"] for align_obj in align_list 
                            if align_obj["has_SA"]])
            mapq = max_mapq
            list_taxid_target = []
            for align_not_suppl in no_suppl_list:
                assert(align_not_suppl["mapq"] == 0)
                align_not_suppl["mapq"] = max_mapq

        list_taxid_target, de_list = [], []
        for alignment in no_suppl_list:
            taxid_target = conv_seqid2taxid[alignment["ref_name"]]
            assert(taxid_target) # If taxid not 'None'
            list_taxid_target.append(taxid_target)
            de_list.append(alignment["de"])
    
        set_taxid_target = set(list_taxid_target)
        if len(set_taxid_target) == 1: # Mergeable directly
            type_alignment = 'second_uniq' 
            taxo_to_write = ';' + str(next(iter(set_taxid_target))) + ';'
        else:
            type_alignment = 'second_plural' 
            taxo_to_write = 's'.join(str(taxid) for taxid in list_taxid_target)
            # print([shared.taxfoo.get_taxid_name(taxid) for taxid in set(list_taxid_target)])

        nb_trashes = sum(map(is_trash, list_taxid_target))
        # print([dico["de"] for dico in no_suppl_list], 
        #       [dico["AS"] for dico in no_suppl_list], 
        #       [dico["mapq"] for dico in no_suppl_list])
        de = max(de_list)
        # de = sum(de_list) / len(de_list)

    else: # Normal case
        type_alignment = 'normal' # i.e. not secondary
        de = representative["de"]
        # if ratio_len >= cutoff_ratio:
        current_taxid = conv_seqid2taxid[representative["ref_name"]]
        nb_trashes = sum(map(is_trash, [current_taxid]))
        taxo_to_write = ';' + str(current_taxid) + ';' # Will be considered as an str


    
    common_to_return += [de]
    # if type_alignment == 'normal':
    #     return [readID, type_alignment, taxo_to_write, ""] + common_to_return
    return [readID, type_alignment, taxo_to_write, nb_trashes] + common_to_return


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
