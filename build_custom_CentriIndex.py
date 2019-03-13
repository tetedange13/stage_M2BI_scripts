#!/usr/bin/env python3

"""

Build a Centrifuge index from a custom database (in the FASTA format

Usage:
  build_custom_CentriIndex.py (-i <inputFaFile>) [-o <dirOut>] [-a <acc2taxid>] [-t <threads>]
  
Options:
  -h --help                help
  --version                version of the script
  -i --infileDB=infile_DB  input file (fa or fq) of the custom database
  -o --dirOut=dir_out      folder where Centrifuge index should be written [default: ./]
  -a --accList=acc_list    file containing the list of the gi in the db [default: none]
  -t --threads=nb_threads  number of threads to use [default: 10]
  
Remarks: 
1) The "-a" parameter becomes mandatory if no file called 'seqid2taxid' 
(mapping sequences IDs with taxids) could be found in 'dirOut/'
2) This file should be formated like that: 'seqid\tGI'
"""

import sys, os, gzip, re
import os.path as osp
import subprocess as sub
import time as t
import multiprocessing as mp
import itertools as ittls
from Bio import Entrez
from functools import partial
from docopt import docopt
import src.remote
import src.check_args as check
import src.parallelized as pll
import src.remote as remote


def map_headers_to_taxid(dirOutput, in_fa_file, taxa_id):
    """
    Very very basic function to map headers against taxid 
    (input file needed by centrifuge-build)
    
    More advanced function can be found in the test.sh file
    """
    with open(seqid2taxid_file, 'a') as outfile:
        for line in in_fa_file:
            if line.startswith(">"):
                outfile.write(line.rstrip('\n').lstrip('>') + '\t' + taxa_id + 
                              '\n')
            

# MAIN:
if __name__ == "__main__":
    # Get arguments:
    ARGS = docopt(__doc__, version='0.1')
    # input_db_path = check.infile(ARGS["--infileDB"], 
    #                              ('fa', 'fasta', 'fq', 'fastq'))[0]
    dirOut = ARGS["--dirOut"]                            
    list_gi_path = ARGS["--accList"]
    NB_THREADS = int(check.input_nb(ARGS["--threads"]))

    # Common variables:
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    #dirIn_db = "/projets/metage_ONT_2019/databases/Zymo_genomes-ZR160406/"
    to_seqid2taxid = osp.join(dirOut, "seqid2taxid")
    conv_headers_Zymo = {
    "Salmonella_enterica_complete_genome":"Salmonella enterica",
    "CP015447.2":"Staphylococcus aureus",
    "Pseudomonas_aeruginosa_complete_genome":"Pseudomonas aeruginosa",
    "Lactobacillus_fermentum_complete_genome":"Lactobacillus fermentum",
    "CP017251.1":"Escherichia coli",
    "LM1":"Listeria monocytogenes",
    "CP015998.1":"Enterococcus faecalis",
    "BS.pilon.polished.v3.ST170922":"Bacillus subtilis",
    "SC":"Saccharomyces cerevisiae",
    "CN":"Cryptococcus neoformans"
    }


    if not osp.isfile(to_seqid2taxid) and list_gi_path == "none":
        print("ERROR! As 'dirOut/seqid2taxid' could not be found, you " +
              "need to specify a file containing the gi of the species " + 
              "contained in your database")
        sys.exit(2)
    
    
    if not osp.isfile(to_seqid2taxid): # GETTING TAXIDS FROM NAME OR GIs:
        id_type = "gi"
        
        if id_type == "sp_name":
            print("Mapping taxids to seqids...")
            with open(list_gi_path, 'r') as headers_file:
                list_headers = headers_file.read().splitlines()

            seqid2taxid_file = open(to_seqid2taxid, 'w')
            dict_seqid2taxid = {}
            for header in list_headers:
                start_header = header[0:2]
                if start_header == "SC" or start_header == "CN":
                # Saccharomyces_cerevisiae or Cryptococcus_neoformans
                    if start_header not in dict_seqid2taxid.keys():
                        sp_name = conv_headers_Zymo[start_header]
                        queried_taxid = remote.query_taxid(sp_name)
                        dict_seqid2taxid[start_header] = queried_taxid

                    seqid2taxid_file.write(header + '\t' + 
                                           dict_seqid2taxid[start_header] + 
                                           '\n')

                else:
                    splitted_header = header.split()[0]
                    sp_name = conv_headers_Zymo[splitted_header]
                    taxid = remote.query_taxid(sp_name)
                    seqid2taxid_file.write(splitted_header + '\t' + taxid + 
                                           '\n')

            seqid2taxid_file.close()


        elif id_type == "gi":
            TAXID_START_TIME = t.time()

            # Get set
            with open(list_gi_path, 'r') as list_gi_file:
               set_acc = {line.rstrip('\n').split()[1] for line in list_gi_file}
            assert(set_acc) # Be sure it is NOT empty

            nb_to_find, nb_found = len(set_acc), 0
            to_acc2taxid = osp.join(dirOut, "acc2taxid")
            to_need_remote = osp.join(dirOut, "need_remote")
            
            if not osp.isfile(to_need_remote) and not osp.isfile(to_acc2taxid):
                # Parallel local taxid mapping:
                to_gz_acc2taxid = (to_dbs + "nt_db/taxo_18feb19/" +
                                   "nucl_wgs.accession2taxid.gz")
                size_chunks = 1000000
                results = []
                pool = mp.Pool(NB_THREADS)
                with gzip.open(to_gz_acc2taxid, 'rt') as acc2taxid_gz:
                    acc2taxid_gz.readline() # Skip header
                    worker = partial(pll.taxid_mapping, set_accession=set_acc)
                    print("\nMapping taxids against accession numbers...")

                    compt_cut_done = 1
                    while True:
                        print("cut", compt_cut_done, ":", 
                              NB_THREADS*size_chunks)
                        compt_cut_done += 1
                        groups = [list(ittls.islice(acc2taxid_gz, size_chunks)) for i in range(NB_THREADS)]

                        # As soon as 1st element not an empty list
                        # (if the 1st one IS empty, the following will beb too)
                        if groups[0]: 
                            result = pool.imap_unordered(worker, groups)
                            results.extend(result)
                        else:
                            print("TOUT CONSOMME")
                            break

                pool.close()
                print("PARALLEL SEARCH TIME:", t.time() - TAXID_START_TIME)

                dict_acc2taxid = {} # To convert accession number to taxid
                set_to_remote = set() # Not found locally, will need remote search
                acc2taxid_file = open(to_acc2taxid, 'w')

                for res in results:
                    if res: # If  NOT 'None'
                        for found_acc_nb, found_taxid in res:
                            dict_acc2taxid[found_acc_nb] = found_taxid
                            acc2taxid_file.write(found_acc_nb + '\t' + 
                                                 found_taxid + '\n')
                acc2taxid_file.close()

                found_acc_numbers = dict_acc2taxid.keys()
                nb_found = len(found_acc_numbers)
                print("TO_FIND:", nb_to_find)
                print("--> FOUND:", nb_found)
                print(len(set_to_remote), "acc_nb couldn't be found locally!")
                print()

                if nb_found == 0:
                    os.remove(to_acc2taxid)

                with open(to_need_remote, 'w') as need_remote_file:
                    for acc_nb in set_acc:
                        if acc_nb not in found_acc_numbers:
                            need_remote_file.write(acc_nb + '\n')
                sys.exit()

            # We solve the problems by querying remotely the taxids
            if osp.isfile(to_need_remote):
                if osp.isfile(to_acc2taxid):
                    dict_acc2taxid = {}
                    with open(to_acc2taxid, 'r') as acc2taxid_file:
                        for line in acc2taxid_file:
                            acc_nb, taxid = line.rstrip('\n').split('\t')
                            dict_acc2taxid[acc_nb] = taxid
                else:
                    dict_acc2taxid = {}

                # Some need to be added manually:
                dict_acc2taxid["CP006690"] = "1051650"
                dict_acc2taxid["CDRP01002148"] = "562"
                dict_acc2taxid["FOMU01000028"] = "321267"
                dict_acc2taxid["CP001098"] = "373903"
                dict_acc2taxid["LIGM01000002"] = "337330"

                set_to_remote = set()
                with open(to_need_remote, 'r') as need_remote_file:
                    for line in need_remote_file:
                        currrent_acc_nb = line.rstrip('\n')
                        if currrent_acc_nb not in dict_acc2taxid.keys():
                            set_to_remote.add(currrent_acc_nb)

                print("REMAINING ONE QUERIED REMOTELY...")
                for pb_acc_nb in set_to_remote:
                    # print("QUERYING REMOTELY:", problematic_acc_nb)
                    queried_taxid = remote.taxid_from_gb_entry(pb_acc_nb)
                    dict_acc2taxid[pb_acc_nb] = queried_taxid
                    with open(to_acc2taxid, 'a') as acc2taxid_file:
                        acc2taxid_file.write(pb_acc_nb + '\t' + 
                                             queried_taxid + '\n')
                
                print("REMOTE SEARCH FINISHED !\n")
                sys.exit()

            
            # Write seqid2taxid file:
            to_seqid2taxid = osp.join(dirOut, "seqid2taxid")
            found_acc_numbers = dict_acc2taxid.keys()

            with open(list_gi_path, 'r') as list_gi_file, \
                 open(to_seqid2taxid, 'w') as seqid2taxid_file:
                for line in list_gi_file:
                    seqid, acc_nb = line.rstrip('\n').split('\t')

                    if acc_nb in found_acc_numbers:
                        seqid2taxid_file.write(seqid + '\t' + 
                                               dict_acc2taxid[acc_nb] + '\n')
                    else:
                        print("STILL A PB..")
                        os.remove(to_seqid2taxid)
                        sys.exit(2)
            

            # os.remove(to_need_remote)
            sys.exit()
                
            print("TOTAL TAXID SEARCH RUNTIME:", t.time() - TAXID_START_TIME)
            
            # Get list of taxids to give to 'centrifuge-download':
            # set_taxids = dict_acc2taxid.values()


    else:
        print("--> FOUND seqid2taxid file")

    # Get list of taxids to give to 'centrifuge-download':
    with open(to_seqid2taxid, 'r') as seqid2taxid_file:
        set_taxids = {line.rstrip('\n').split('\t')[1] for line in seqid2taxid_file}


    # GET TAXONOMY FILES:
    if not (osp.isfile(dirOut+"nodes.dmp") and osp.isfile(dirOut+"names.dmp")):
        print("\nGenerating taxonomic tree and names...")
        with open("toto.txt", 'w') as toto:
            toto.write(','.join(list(set_taxids)))
        
        cmd_centriDl = ("centrifuge-download -v -o " + dirOut + " -t " + 
                        ','.join(list(set_taxids)) + " taxonomy")
        sub.Popen(cmd_centriDl.split(), shell=True)
        print("Done !\n")
    
    else:
        print("--> FOUND all taxonomic files in", dirOut, "!\n")
    
    
    # BUILDING INDEX
    path_to_centriBuild = ("/home/sheldon/Applications/centrifuge-master21" +
                               "jan19/build/bin/centrifuge-build")
    db_name = "silva"
    cmd_centriBuild = (path_to_centriBuild + " -p 7 --conversion-table " + 
                       dirOut + "seqid2taxid " + "--taxonomy-tree " + 
                       dirOut + "nodes.dmp --name-table " + dirOut + 
                       "names.dmp " + input_db_path + " " + dirOut + db_name)
    print("Building centrifuge index...")
    with open(dirOut + "centri-build.log", 'w') as centriBuild_log:
        sub.Popen(cmd_centriBuild.split(), stdout=centriBuild_log).communicate()
        
    print("Index building finished !")
    print("A log file can be found in", dirOut, '\n')
    