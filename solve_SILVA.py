#!/usr/bin/env python3


import sys, os, gzip
import os.path as osp
import time as t
import multiprocessing as mp
import itertools as ittls
from functools import partial
import src.remote as remote
import src.parallelized as pll

# import src.shared_fun as shared


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


def test(to_wrong2good_file, to_seqid2taxid):
    with open(to_wrong2good_file, 'r') as wrong2good_file:
        set_found_taxids = {line.split('\t')[0] for line in wrong2good_file}

    done_taxids = set()
    with open(to_seqid2taxid, 'r') as seqid2taxid_file, \
         open("problems.txt", 'w') as pbs_file:
        for line in seqid2taxid_file:
            seqid, taxid = line.rstrip('\n').split('\t')
            if taxid not in set_found_taxids and taxid not in done_taxids:
                pbs_file.write(taxid + '\t' + seqid.split('.')[0] + '\n')
                done_taxids.add(taxid)



def detect_problems(to_seqid2taxid, to_seqid2acc, to_target_list, dirOut="./"):
    """
    to_target_list = list of seqid that have been mapped
    """
    import src.ncbi_taxdump_utils as taxo_utils
    print("Loading taxonomic module...")
    taxfoo = taxo_utils.NCBI_TaxonomyFoo()
    nodes_path = to_dbs + "nt_db/taxo_18feb19/nodes.dmp"
    names_path = to_dbs + "nt_db/taxo_18feb19/names.dmp"
    taxfoo.load_nodes_dmp(nodes_path)
    taxfoo.load_names_dmp(names_path)
    print("Module loaded !")
    print(taxfoo.get_lineage(161659));sys.exit()

    dict_seqid2taxid = {}
    with open(to_seqid2taxid, 'r') as seqid2taxid_file:#, \
         #open("problems.txt", 'w') as pbs_file:
        for line in seqid2taxid_file:
            seqid, taxid = line.rstrip('\n').split('\t')
            assert(seqid not in dict_seqid2taxid.keys())
            dict_seqid2taxid[seqid] = taxid
            # print(list(dict_seqid2taxid.keys())[list(dict_seqid2taxid.values()).index(acc_nb)])    

    dict_seqid2acc = {}
    with open(to_seqid2acc) as seqid2acc_file:
        for line in seqid2acc_file:
            seqid, acc_nb = line.rstrip('\n').split()
            dict_seqid2acc[seqid] = acc_nb


    set_problems = set()
    done_taxids = set()
    set_uniden = set()
    with open(to_target_list, 'r') as target_list:
        for line in target_list:
            ref_name = line.rstrip('\n')
            current_taxid = dict_seqid2taxid[ref_name]

            name_taxid = taxfoo.get_taxid_name(int(current_taxid))
            if name_taxid == "unidentified":
                set_uniden.add(current_taxid)
                # sys.exit()
            if not name_taxid and current_taxid not in done_taxids:
            # if not name_taxid:
                set_problems.add(current_taxid + '\t' + 
                                 dict_seqid2acc[ref_name] + '\n')
                done_taxids.add(current_taxid)

    if set_problems:
        print("PROBLEMS FOUND..")
        with open(osp.join(dirOut, "problems.txt"), 'w') as pbs_file:
            [pbs_file.write(problem) for problem in set_problems]
    else:
        print("NO PROBLEMS FOUND !")
    print("Taxids assiciated to 'unidentified':", set_uniden)
    print()


def remote_taxid_search(to_problems_file, dirOut="./"):
    """
    """
    print("REMOTE TAXID SEARCH...")
    dict_wrong2good = {} # Conversion between wrong and good taxid
    to_wrong2good_file = osp.join(dirOut, "wrong2good_taxids")

    if osp.isfile(to_wrong2good_file):
        with open(to_wrong2good_file, 'r') as wrong2good_file:
            for line in wrong2good_file:
                wrong_taxid, good_taxid = line.rstrip('\n').split('\t')
                dict_wrong2good[wrong_taxid] = good_taxid 

    possible_issues = ("HTTP400", "PB_TAXID_SEARCH", "NOT_FOUND")
    with open(to_problems_file, 'r') as pbs_file:
        for line in pbs_file:
            wrong_taxid, acc_nb = line.rstrip('\n').split('\t')

            if wrong_taxid not in dict_wrong2good.keys():
                good_taxid = remote.taxid_from_gb_entry(acc_nb)

                if good_taxid in possible_issues:
                    with open(osp.join(dirOut, "remote_fails"), 'w') as fails:
                        fails.write(acc_nb + '\t' + good_taxid + '\n')
                else:
                    dict_wrong2good[wrong_taxid] = good_taxid
                    with open(to_wrong2good_file, 'a') as wrong2good_file:
                        wrong2good_file.write(wrong_taxid + '\t' + 
                                      dict_wrong2good[wrong_taxid] + '\n')
    print("FINISHED !\n")


def local_taxid_search(gz_db_name, to_problems_file, dirOut="./", 
                       nb_threads=10):
    """
    Parallel local taxid mapping:
    """
    to_dbs = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
              "metage_ONT_2019/nt_db/taxo_18feb19/")
    to_gz_acc2taxid = to_dbs + "nucl_" + gz_db_name + ".accession2taxid.gz"
    TAXID_START_TIME = t.time()
    dict_acc2wrong = {}

    with open(to_problems_file, 'r') as pbs_file:
        for line in pbs_file:
            wrong_taxid, acc_nb = line.rstrip('\n').split('\t')
            if wrong_taxid not in dict_acc2wrong.values():
                dict_acc2wrong[acc_nb] = wrong_taxid

    size_chunks = 1000000
    results = []
    my_pool = mp.Pool(nb_threads)
    with gzip.open(to_gz_acc2taxid, 'rt') as acc2taxid_gz:
        acc2taxid_gz.readline() # Skip header
        worker = partial(pll.taxid_mapping, 
                         set_accession=set(dict_acc2wrong.keys()))
        print("\nMapping taxids against accession numbers...")
        print("DB USED:", gz_db_name)

        compt_cut_done = 1
        while True:
            print("cut", compt_cut_done, ":", nb_threads*size_chunks)
            compt_cut_done += 1
            groups = [list(ittls.islice(acc2taxid_gz, size_chunks)) for i in range(nb_threads)]

            # As soon as 1st element not an empty list
            # (if the 1st one IS empty, the following will beb too)
            if groups[0]: 
                result = my_pool.imap_unordered(worker, groups)
                results.extend(result)
            else:
                print("TOUT CONSOMME")
                break
    my_pool.close()

    dict_acc2good = {} # To convert accession number to GOOD taxid
    to_wrong2good_file = os.path.join(dirOut, "wrong2good_taxids")
    with open(to_wrong2good_file, 'a') as wrong2good_file:
        for res in results:
            if res: # If  NOT 'None'
                for found_acc_nb, found_taxid in res:
                    if found_taxid not in dict_acc2good.values():
                        old_taxid = dict_acc2wrong[found_acc_nb]
                        dict_acc2good[found_acc_nb] = found_taxid
                        wrong2good_file.write(old_taxid + '\t' + 
                                              found_taxid + '\n')

    nb_found = len(dict_acc2good.keys())
    nb_to_find = len(dict_acc2wrong.keys())
    if nb_found != 0:
        to_need_remote = os.path.join(dirOut, gz_db_name +"_need_remote")
        with open(to_need_remote, 'w') as need_remote_file:
            found_acc_numbers = dict_acc2good.keys()
            for acc_to_find in dict_acc2wrong:
                if acc_to_find not in found_acc_numbers:
                    need_remote_file.write(dict_acc2wrong[acc_to_find] + '\t' + 
                                           acc_to_find + '\n')

    print("PARALLEL SEARCH TIME:", t.time() - TAXID_START_TIME)
    print("TO_FIND:", nb_to_find)
    print("--> FOUND:", nb_found)
    print(nb_to_find - nb_found, "acc_nb couldn't be found locally!")
    print()


def correct_seqid2taxid(to_old_seqid2taxid, to_wrong2good_taxids, dirOut='./'):
    """
    Correct seqid2taxid file, by replacing taxids that have been found to be
    wrong (taxid not found in the nt NCBI taxonomy dump files), by correct 
    taxids (possibly found remotely or locally)
    """
    print("CORRECTING TAXIDS...")
    dict_wrong2good = {}
    with open(to_wrong2good_taxids, 'r') as wrong2good_file:
        for line in wrong2good_file:
            wrong_taxid, good_taxid = line.rstrip('\n').split('\t')
            assert(wrong_taxid not in dict_wrong2good.keys())
            dict_wrong2good[wrong_taxid] = good_taxid

    all_wrong_taxids = dict_wrong2good.keys()
    compt_corr = []
    with open(to_old_seqid2taxid, 'r') as old_seqid2taxid, \
         open(osp.join(dirOut, "seqid2taxid"), 'w') as new_seqid2taxid:
        for line in old_seqid2taxid:
            seqid, old_taxid = line.rstrip('\n').split('\t')
            if old_taxid in all_wrong_taxids:
                compt_corr.append(old_taxid)
                new_seqid2taxid.write(seqid + '\t' + 
                                       dict_wrong2good[old_taxid] + '\n')
            else:
                new_seqid2taxid.write(line)

    print("FINISHED!", len(compt_corr), "CORRECTIONS DONE")
    print("(representing", len(set(compt_corr)), "unique taxids)")
    print()



# MAIN:
if __name__ == "__main__":
    to_dbs = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
              "metage_ONT_2019/")
    to_dbs_SILVA = to_dbs + "SILVA_refNR132_28feb19/"

    # test(to_dbs_SILVA + "wrong2good_taxids", "seqid2taxid")

    # detect_problems("seqid2taxid", to_dbs_SILVA + "seqid2acc", 
    #                 "bonj/headers.txt")

    # local_taxid_search("wgs", "problems.txt")
    # local_taxid_search("gb", "wgs_need_remote")
    # local_taxid_search("gss", "gb_need_remote")

    remote_taxid_search("gb_need_remote")
    
    # correct_seqid2taxid(to_dbs + "Centri_idxes/SILVA/seqid2taxid", 
    #                     "bonj/wrong2good_taxids")

    #solve_headers(to_dbs_SILVA + "taxmap_ncbi_ssu_ref_nr99_132.txt", 
    #              to_dbs_SILVA + "headers.txt")
