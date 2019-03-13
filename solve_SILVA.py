#!/usr/bin/env python3


import sys, os
import os.path as osp
import src.remote as remote
import src.ncbi_taxdump_utils as taxo_utils


def solve_headers(to_taxmap_file, to_headers_list_file):
    """
    """
    dict_acc2taxid = {}
    with open(to_taxmap_file, 'r') as taxmap:
        for line in taxmap:
            acc_nb, _, taxid = line.rstrip('\n').split('\t')[0:3]
            
            if acc_nb not in dict_acc2taxid.keys():
                dict_acc2taxid[acc_nb] = taxid
            else:
                pass
                #print("DUPLICATE:", acc_nb, taxid)
                #sys.exit(2)
        
    my_set = set()
    with open(to_headers_list_file, 'r') as headers_file, \
         open("seqid2taxid", 'w') as seqid2taxid_file:
        for header in headers_file:
            splitted_header = header.rstrip('\n').split('.')[0] # We get only the 1st part
            if header not in my_set:
                my_set.add(header)
                seqid2taxid_file.write(header.rstrip('\n') + '\t' + 
                                       dict_acc2taxid[splitted_header] + '\n')
            else:
                print("PROBLEM")
                print(header)
                sys.exit(2)
      

# def eval_align(ref_name, dict_conv_seqid2taxid, taxonomic_cutoff):
def eval_align(current_taxid, ref_name, set_problems):
    """
    """
    # current_taxid = dict_conv_seqid2taxid[ref_name]
    sp_name = taxfoo.get_taxid_name(current_taxid)
    set_problems = set()

    if not sp_name:
        # .add(ref_name, current_taxid)
        print(ref_name, "HAS A PROBLEM !");sys.exit()

        # to_problems_file = "problematic_ref_names.txt"
        # if osp.isfile(to_problems_file):
        #   with open(to_problems_file, 'r') as problems:
        #     set_taxids = {line.rstrip('\n').split('\t')[1] for line in problems}
        # else:
        #   set_taxids = set()

        # if str(current_taxid) not in set_taxids:
        #     queried_taxid = remote.query_taxid(ref_name) # Try remotely
        #     if queried_taxid == "TAXO_NOT_FOUND": # Still not working
        #         with open(to_problems_file, 'a') as problems:
        #             problems.write(ref_name + '\t' + str(current_taxid) + '\n')
        #         print("PROBLEM with:", ref_name, "of associated taxid:", current_taxid)
        #         return None
        #     else:
        #         print("LE REMOTE A MARCHE !")
        #         sp_name = taxfoo.get_taxid_name(queried_taxid)

        # else:
        #   return None

    # assert(sp_name) # Not 'None'
    #current_rank = taxfoo.get_taxid_rank(current_taxid)
    #assert(current_rank) # Not 'None'

    #return in_zymo(sp_name, current_rank, taxonomic_cutoff)


def solve_pb_taxids(to_old_seqid2taxid):
    """
    """
    set_problems = set()
    with open(to_old_seqid2taxid, 'r') as old_seqid2taxid, \
         open("problems.txt", 'w') as pbs_file:
        for line in old_seqid2taxid:
            seqid, old_taxid = line.rstrip('\n').split('\t')

            if not taxfoo.get_taxid_name(int(old_taxid)):
                print(seqid, old_taxid, "HAS A PROBLEM !")

                if old_taxid not in set_problems:
                    pbs_file.write(seqid.split('.')[0] + ' ' + old_taxid + '\n')
                    set_problems.add(old_taxid)
                    # new_taxid = remote.taxid_from_gb_entry()
                    # print(new_taxid);sys.exit()
            # res_eval = eval_align(old_taxid, seqid)
            # if res_eval:
            #     list_evals.append(res_eval)
            # dict_seqid2taxid[seqid] = old_taxid

    print(len(set_problems))
    sys.exit()

    to_toto = "tempo.txt"
    wrong2good_taxid = {}
    
    if osp.isfile(to_toto):
        with open(to_toto) as toto:
             for line in toto:
                false_taxid, good_taxid = line.rstrip('\n').split('\t')
                wrong2good_taxid[false_taxid] = good_taxid
    
    if True:
        with open("problematic_ref_names.txt", 'r') as pb_ref_names_file:
            for line in pb_ref_names_file:
                ref_name_to_split, wrong_taxid = line.rstrip('\n').split('\t')
                
                if wrong_taxid not in wrong2good_taxid.keys():
                    ref_name = ref_name_to_split.split('.')[0]
                    print("REQUESTING:", ref_name, "..")
                    if ref_name == "MCIE01000002":
                        taxid = 46170
                    elif ref_name == "JTBM01000001":
                        taxid = 28901
                    else:
                        taxid = remote.taxid_from_gb_entry(ref_name)
                        wrong2good_taxid[wrong_taxid] = str(taxid)
                    
                    with open(to_toto, 'a') as toto:
                        toto.write(wrong_taxid + '\t' + str(taxid) + '\n')

    
    # Correct seqid2taxid file:
    all_wrong_taxids = wrong2good_taxid.keys()
    with open(to_old_seqid2taxid, 'r') as old_seqid2taxid, \
         open("seqid2taxid", 'w') as seqid2taxid_file:
        for line in old_seqid2taxid:
            seqid, old_taxid = line.rstrip('\n').split('\t')
            if old_taxid in all_wrong_taxids:
                seqid2taxid_file.write(seqid + '\t' + 
                                       wrong2good_taxid[old_taxid] + '\n')
            else:
                seqid2taxid_file.write(line) 
    

# MAIN:
if __name__ == "__main__":
    to_dbs = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
              "metage_ONT_2019/")
    to_dbs_SILVA = to_dbs + "SILVA_refNR132_28feb19/"
    
    #solve_headers(to_dbs_SILVA + "taxmap_ncbi_ssu_ref_nr99_132.txt", 
    #              to_dbs_SILVA + "headers.txt")

    global taxfoo
    taxfoo = taxo_utils.NCBI_TaxonomyFoo()
    nodes_path = to_dbs + "nt_db/taxo_18feb19/nodes.dmp"
    names_path = to_dbs + "nt_db/taxo_18feb19/names.dmp"
    taxfoo.load_nodes_dmp(nodes_path)
    taxfoo.load_names_dmp(names_path)
    solve_pb_taxids(to_dbs + "Centri_idxes/SILVA/seqid2taxid")
    
