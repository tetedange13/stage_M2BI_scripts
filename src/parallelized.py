#!/usr/bin/env python3


"""
Contains all functions used to parallelise processes
"""


from sys import exit
import src.taxo_eval as eval
Series = eval.pd.Series


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
            
        if 'AS' in representative.keys():
            max_AS = max(map(lambda dico: dico["AS"], no_suppl_list))
            only_equiv_list = [dico for dico in no_suppl_list 
                               if dico["AS"] == max_AS]

            if representative["has_SA"]:
                # print(readID, "suppl_as_repr")
                max_mapq = max([align_obj["mapq"] for align_obj in align_list 
                                if align_obj["has_SA"]])
                returned_dict['mapq'] = max_mapq
                # for align_not_suppl in no_suppl_list:
                for align_not_suppl in only_equiv_list:
                    # Need to propagate max value of MAPQ to all secondaries:
                    align_not_suppl["mapq"] = max_mapq
                del align_not_suppl

        else:
            only_equiv_list = no_suppl_list
        
        list_taxid_target = [conv_seqid2taxid[align_dict["ref_name"]]
                             for align_dict in only_equiv_list]
        assert(all(list_taxid_target)) # Check if there are NO 'None'

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
    taxo_to_write = (';' + 's'.join(str(taxid) for taxid 
                                               # in sorted(list_taxid_target)) +
                                               in list_taxid_target) +
                     ';')
    # return [readID, type_alignment, taxo_to_write, nb_trashes, mapq, len_align,
    #         ratio_len, de]
    returned_dict['lineage'] = taxo_to_write
    returned_dict['nb_trashes'] = nb_trashes
    returned_dict['type_align'] = type_alignment
    return returned_dict


# def eval_taxo(one_csv_index_val, two_col_from_csv, set_levels_prok, 
#               taxonomic_cutoff, mode):
def eval_taxo(one_line_from_csv, set_levels_prok, taxonomic_cutoff, mode):
    """
    Do the parallel taxonomic evaluation on a given Pandas DataFrame
    The 'mode' parameter allows to change between LCA and majo voting, for the
    handling of multi-hits/mult-mapping
    """
    # /!\ CAREFUL WITH THE ORDER HERE:
    one_csv_index_val, lineage_val, type_align = one_line_from_csv.values
    remark_eval = type_align

    if type_align == 'normal':
        taxid_to_eval = int(lineage_val.strip(';'))
        
    elif 'second' in type_align: # Secondary (with 1 unique taxid or more)
        list_taxid_target = list(map(int, lineage_val.strip(';').split('s')))
        if type_align == 'second_uniq':
            taxid_to_eval = list_taxid_target[0]

        else: # More than 1 unique taxid
            if mode == 'MAJO_OLD':
                res_second_handling = eval.majo_voting_old(list_taxid_target, 'species')
            elif mode == 'MINOR_RM_LCA':
                res_second_handling = eval.majo_voting(list_taxid_target)
            elif mode == 'LCA':
                res_second_handling = eval.make_lca(list_taxid_target)
            elif mode == 'TOP_ONE':
                # Like that, we always take the 1st one
                res_second_handling = ('topOne_only', list_taxid_target[0])
            else:
                print("ERROR: WRONG MODE !")
                sys.exit()
            remark_eval = res_second_handling[0]

            if False:
            # if remark_eval == 'no_majo_found': # Problem (
                remark_eval += ';lca'
                taxid_to_eval = eval.taxfoo.find_lca(set(list_taxid_target))
            elif len(res_second_handling) == 1: # Other problem ('no_majo_found' mainly)
                remark_eval = type_align + ';' + remark_eval
                # /!\ CAREFUL WITH THE ORDER HERE:
                return Series([one_csv_index_val, 'no_majo_found', 
                                       'no_majo_found', 'FP', remark_eval],
                                      index=['index', 'taxid_ancester', 
                                             'final_taxid', 'res_eval', 
                                             'remark_eval'])
            else:
                taxid_to_eval = res_second_handling[1]

    else:
        print("ERROR ! 'type_align' NOT KNOWN !")
        exit()
    #     taxid_to_eval = lineage_val # Should be a nan
    #     assert(pd.np.isnan(taxid_to_eval))

    taxid_ancester, classif, remark = eval.in_zymo(taxid_to_eval, 
                                                   set_levels_prok, 
                                                   taxonomic_cutoff)
    remark_eval += ';' + remark

    # /!\ CAREFUL WITH THE ORDER HERE:
    return Series([one_csv_index_val, taxid_ancester, int(taxid_to_eval), 
                   classif, remark_eval], 
                  index=['index', 'taxid_ancester', 'final_taxid', 'res_eval', 
                  'remark_eval'])



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
