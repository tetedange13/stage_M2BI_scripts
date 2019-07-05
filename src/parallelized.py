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
    Determine wheteher a given taxid is a 'trash', i.e. a poorly characterized
    entry (mainly for SILVA) 
    """
    taxid_name = eval.taxfoo.get_taxid_name(int(taxid_to_eval))
    if not taxid_name:
        print('WARNING: {} is an unknown taxid !'.format(taxid_to_eval))
        return False   
    
    tax_name_to_eval = taxid_name.lower()

    if 'metagenome' in tax_name_to_eval: return True
    elif 'uncultured' in tax_name_to_eval:
        return True 
    elif 'unidentified' in tax_name_to_eval: return True
    elif 'phage' in tax_name_to_eval: return True
    elif tax_name_to_eval == 'synthetic construct': return True
    elif 'virus' in tax_name_to_eval: return True
        
    return False


def eval_taxo(one_line_from_csv, set_levels_prok, taxonomic_cutoff, mode):
    """
    Do the parallel taxonomic evaluation on a given Pandas DataFrame
    The 'mode' parameter allows to change between LCA and majo voting, for the
    handling of multi-hits/mult-mapping
    """
    # /!\ CAREFUL WITH THE ORDER HERE:
    one_csv_index_val, lineage_val, type_align = one_line_from_csv.values
    remark_eval = type_align

    if len(lineage_val) == 1:
        taxid_to_eval = lineage_val[0]
        
    else: # Secondary (with 1 unique taxid or more)
        list_taxid_target = lineage_val
        if type_align == 'second_uniq':
            taxid_to_eval = list_taxid_target[0]

        else: # More than 1 unique taxid
            if mode == 'MAJO_OLD':
                res_second_handling = eval.majo_voting_old(list_taxid_target, 
                                                           'species')
            elif mode == 'MINOR_RM_LCA':
                res_second_handling = eval.majo_voting(list_taxid_target)
            elif mode == 'LCA':
                res_second_handling = eval.make_lca(set(list_taxid_target))
            elif mode == 'TOP_ONE':
                # Like that, we always take the 1st one
                res_second_handling = ('topOne_only', list_taxid_target[0])
            else:
                print("ERROR: WRONG MODE !")
                sys.exit()
            remark_eval = res_second_handling[0]

            if len(res_second_handling) == 1: # Other problem ('no_majo_found' mainly)
                remark_eval = type_align + ';' + remark_eval
                # /!\ CAREFUL WITH THE ORDER HERE:
                return Series([one_csv_index_val, 'no_majo_found', 
                               'no_majo_found', 'miss', remark_eval],
                              index=['readID', 'taxid_ancester', 'final_taxid', 
                                     'res_eval', 'remark_eval'])
            else:
                taxid_to_eval = res_second_handling[1]

    taxid_ancester, classif, remark = eval.in_zymo(taxid_to_eval, 
                                                   set_levels_prok, 
                                                   taxonomic_cutoff)
    remark_eval += ';' + remark

    # /!\ CAREFUL WITH THE ORDER HERE:
    return Series([one_csv_index_val, taxid_ancester, int(taxid_to_eval), 
                   classif, remark_eval], 
                  index=['readID', 'taxid_ancester', 'final_taxid', 'res_eval', 
                         'remark_eval'])



# FOR 0-SOLVE_SILVA.PY:
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
