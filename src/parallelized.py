#!/usr/bin/env python3


"""
Contains all functions used to parallelise processes
"""


import src.taxo_eval as eval


# FROM STAT_TAX.PY:
def is_trash(taxid_to_eval):
    """
    """
    taxid_name = eval.taxfoo.get_taxid_name(int(taxid_to_eval))
    if not taxid_name:
        print(taxid_to_eval)
    assert(taxid_name)
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


        repr_AS = no_suppl_list[0]["AS"]
        toto = [(idx, dico["AS"]) for idx, dico in enumerate(no_suppl_list) 
                if dico["AS"] == repr_AS and idx != 0]
        # Contains ONLY equivalent (same AS score) align
        only_equiv_list = [dico for dico in no_suppl_list 
                           if dico["AS"] == repr_AS]
        if toto:
            print("LAST_EQUIV:", toto[-1][0])
        else:
            print("NO_EQUIV")
        # Need to propagate max value of MAPQ to all secondaries:
        if representative["has_SA"]:
            # print(readID, "suppl_as_repr")
            max_mapq = max([align_obj["mapq"] for align_obj in align_list 
                            if align_obj["has_SA"]])
            mapq = max_mapq
            list_taxid_target = []
            for align_not_suppl in only_equiv_list:
                assert(align_not_suppl["mapq"] == 0)
                align_not_suppl["mapq"] = max_mapq
            del align_not_suppl

        list_taxid_target, de_list = [], []
        for alignment in only_equiv_list:
            taxid_target = conv_seqid2taxid[alignment["ref_name"]]
            assert(taxid_target) # If taxid not 'None'
            list_taxid_target.append(taxid_target)
            de_list.append(alignment["de"])
        del alignment
    
        if len(set(list_taxid_target)) == 1: # Mergeable directly
            type_alignment = 'second_uniq' 
        else:
            type_alignment = 'second_plural'
            # print([shared.taxfoo.get_taxid_name(taxid) for taxid in set(list_taxid_target)])
        
        taxo_to_write = 's'.join(str(taxid) for taxid in list_taxid_target)          
        nb_trashes = sum(map(is_trash, list_taxid_target))
        # print([dico["de"] for dico in no_suppl_list], 
        #       [dico["AS"] for dico in no_suppl_list], 
        #       [dico["mapq"] for dico in no_suppl_list])
        de = max(de_list)
        # de = sum(de_list) / len(de_list)

    else: # Normal case
        type_alignment = 'normal' # i.e. not secondary
        de = representative["de"]
        current_taxid = conv_seqid2taxid[representative["ref_name"]]
        nb_trashes = int(is_trash(current_taxid))
        taxo_to_write = ';' + str(current_taxid) + ';' # Will be considered as an str

    common_to_return += [de]
    return [readID, type_alignment, taxo_to_write, nb_trashes] + common_to_return


def eval_taxo(one_csv_index_val, two_col_from_csv, sets_levels, 
              taxonomic_cutoff, mode):
    """
    Do the parallel taxonomic evaluation on a given Pandas DataFrame
    The 'mode' parameter allows to change between LCA and majo voting, for the
    handling of multi-hits/mult-mapping
    """
    lineage_val = two_col_from_csv.loc[one_csv_index_val, "lineage"]
    type_align = two_col_from_csv.loc[one_csv_index_val, "type_align"]
    remark_eval = type_align

    # if lineage_val.startswith(';'): # Normal
    if type_align == 'normal':
        taxid_to_eval = lineage_val.strip(';')
        
    else: # Secondary (with 1 unique taxid or more)
        list_taxid_target = list(map(int, lineage_val.split('s')))
        if type_align == 'second_uniq':
            taxid_to_eval = list_taxid_target[0]
        else: # More than 1 unique taxid
            if mode == 'MAJO':
                res_second_handling = eval.majo_voting(list_taxid_target, 
                                                       taxonomic_cutoff)
            elif mode == 'LCA':
                res_second_handling = eval.make_lca(list_taxid_target)
            remark_eval = res_second_handling[0]
            if len(res_second_handling) == 1: # Problem (root reached)
                # list_sp_name = '&'.join([eval.taxfoo.get_taxid_name(taxid) 
                #                          for taxid in list_taxid_target])
                # return (one_csv_index_val, list_sp_name, 'FP', remark_eval)
                return (one_csv_index_val, remark_eval, 'FP', remark_eval)
            else:
                taxid_to_eval = res_second_handling[1]
                if remark_eval == 'majo_notInKey':
                    return (one_csv_index_val, 
                            eval.taxfoo.get_taxid_name(int(taxid_to_eval)), 
                            'FP', remark_eval)

    taxo_name, classif, remark = eval.in_zymo(taxid_to_eval, sets_levels, 
                                              taxonomic_cutoff)
    remark_eval += ';' + remark

    return (one_csv_index_val, taxo_name, classif, remark_eval)



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
    del line

    if to_return:
        return to_return
