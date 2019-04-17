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
        print('WARNING: {} is an unknown taxid !'.format(taxid_to_eval))
        return False   
    # assert(taxid_name)
    
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
    # mapq, ratio_len = representative["mapq"], representative["ratio_len"]
    # len_align = representative["len_align"]
    returned_dict = {'readID' : readID,
                     'mapq' : representative["mapq"], 
                     'ratio_len' : representative["ratio_len"],
                     'len_align' : representative["len_align"]}
    

    if nb_alignments_for_readID > 1: # Secondary alignment
        assert(not representative["is_suppl"])
        assert(not representative["is_second"])

        # Get rid of supplementaries:
        no_suppl_list = [align_obj for align_obj in align_list 
                         if not align_obj["has_SA"]]
        # Check if they are all secondaries (1st one can be the representative):
        # assert(all(map(lambda dico: not dico["is_suppl"], no_suppl_list)))
        # assert(all(map(lambda dico: dico["is_second"], no_suppl_list[1: ])))

        if not no_suppl_list: # Only supplementaries entries for this read
            return {'readID':readID, 'type_align':'only_suppl', 
                    'lineage':'no'}
            # return [readID, 'only_suppl', 'no']

        max_AS = max(map(lambda dico: dico["AS"], no_suppl_list))
        only_equiv_list = [dico for dico in no_suppl_list 
                           if dico["AS"] == max_AS]

        list_taxid_target = []
        if representative["has_SA"]:
            # print(readID, "suppl_as_repr")
            max_mapq = max([align_obj["mapq"] for align_obj in align_list 
                            if align_obj["has_SA"]])
            returned_dict['mapq'] = max_mapq
            # for align_not_suppl in no_suppl_list:
            for align_not_suppl in only_equiv_list:
                assert(align_not_suppl["mapq"] == 0)
                # Need to propagate max value of MAPQ to all secondaries:
                align_not_suppl["mapq"] = max_mapq
                taxid_target = conv_seqid2taxid[align_not_suppl["ref_name"]]
                assert(taxid_target) # If taxid not 'None'
                list_taxid_target.append(taxid_target)
            del align_not_suppl
        else:
            # for align_not_suppl in no_suppl_list:
            for align_not_suppl in only_equiv_list:
                taxid_target = conv_seqid2taxid[align_not_suppl["ref_name"]]
                assert(taxid_target) # If taxid not 'None'
                list_taxid_target.append(taxid_target)
            del align_not_suppl
        
        # if 'de' in align_obj[0].keys():
        #     de_list = [a_dict['de'] for a_dict in align_not_suppl]
        #     returned_dict['de'] = max(de_list)
            # de = sum(de_list) / len(de_list)
        if len(set(list_taxid_target)) == 1: # Mergeable directly
            type_alignment = 'second_uniq' 
        else:
            type_alignment = 'second_plural'
        nb_trashes = sum(map(is_trash, list_taxid_target))

        # print([dico["de"] for dico in no_suppl_list], 
        #       [dico["AS"] for dico in no_suppl_list], 
        #       [dico["mapq"] for dico in no_suppl_list])
        
    else: # Normal case
        type_alignment = 'normal' # i.e. not secondary
        # de = representative["de"]
        current_taxid = conv_seqid2taxid[representative["ref_name"]]
        list_taxid_target = [current_taxid]
        nb_trashes = int(is_trash(current_taxid))
    
     # Taxid has to be considered as an str:
    taxo_to_write = (';' + 
                     's'.join(str(taxid) for taxid in list_taxid_target) + ';')
    # return [readID, type_alignment, taxo_to_write, nb_trashes, mapq, len_align,
    #         ratio_len, de]
    returned_dict['lineage'] = taxo_to_write
    returned_dict['nb_trashes'] = nb_trashes
    returned_dict['type_align'] = type_alignment
    return returned_dict


def eval_taxo(one_csv_index_val, two_col_from_csv, set_levels_prok, 
              taxonomic_cutoff, mode):
    """
    Do the parallel taxonomic evaluation on a given Pandas DataFrame
    The 'mode' parameter allows to change between LCA and majo voting, for the
    handling of multi-hits/mult-mapping
    """
    lineage_val = two_col_from_csv.loc[one_csv_index_val, "lineage"]
    type_align = two_col_from_csv.loc[one_csv_index_val, "type_align"]
    remark_eval = type_align

    if type_align == 'normal':
        taxid_to_eval = lineage_val.strip(';')
        
    else: # Secondary (with 1 unique taxid or more)
        list_taxid_target = list(map(int, lineage_val.strip(';').split('s')))
        if type_align == 'second_uniq':
            taxid_to_eval = list_taxid_target[0]
        else: # More than 1 unique taxid
            if mode == 'MAJO':
                res_second_handling = eval.majo_voting(list_taxid_target, 
                                                       'species')
            elif mode == 'LCA':
                res_second_handling = eval.make_lca(list_taxid_target)
            remark_eval = res_second_handling[0]
            if len(res_second_handling) == 1: # Problem ('no_majo_found')
                # list_sp_name = '&'.join([eval.taxfoo.get_taxid_name(taxid) 
                #                          for taxid in list_taxid_target])
                # return (one_csv_index_val, list_sp_name, 'FP', remark_eval)
                return (one_csv_index_val, 'no_majo_found', remark_eval, 'FP', 
                        remark_eval)
            else:
                taxid_to_eval = res_second_handling[1]
                if remark_eval == 'majo_notInKey':
                    return (one_csv_index_val, int(taxid_to_eval),
                            eval.taxfoo.get_taxid_name(int(taxid_to_eval)), 
                            'FP', remark_eval)

    taxo_name, classif, remark = eval.in_zymo(taxid_to_eval, set_levels_prok, 
                                              taxonomic_cutoff)
    remark_eval += ';' + remark

    return (one_csv_index_val, int(taxid_to_eval), taxo_name, classif, 
            remark_eval)



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
