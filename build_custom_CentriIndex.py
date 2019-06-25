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

import sys, os, re
import os.path as osp
import subprocess as sub
import time as t
from Bio import Entrez
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
    input_db_path = check.infile(ARGS["--infileDB"], 
                                 ('fna', 'fa', 'fasta', 'fq', 'fastq'))[0]
    dirOut = ARGS["--dirOut"]                            
    list_gi_path = ARGS["--accList"]
    NB_THREADS = check.input_nb(ARGS["--threads"]) # TYPE STRING

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
            print("MARCHE PLUS !\n")
            sys.exit()


    else:
        print("--> FOUND seqid2taxid file in", dirOut, "!")

    # Get list of taxids to give to 'centrifuge-download':
    # with open(to_seqid2taxid, 'r') as seqid2taxid_file:
    #     set_taxids = {line.rstrip('\n').split('\t')[1] for line in seqid2taxid_file}


    # GET TAXONOMY FILES:
    if (not (osp.isfile(dirOut+"nodes.dmp") and
        osp.isfile(dirOut + "names.dmp"))):
        print("ERROR: Taxonomic files (.dmp) NOT FOUND in", dirOut, "\n")
        sys.exit()

        print("\nGenerating taxonomic tree and names...")  
        # cmd_centriDl = ("centrifuge-download -v -o " + dirOut + " -t " + 
        #                 ','.join(list(set_taxids)) + " taxonomy")
        # sub.Popen(cmd_centriDl.split(), shell=True)
        cmd_centriDl = "centrifuge-download -v -o ./ taxonomy"
        os.system(cmd_centriDl)
        print("Done !\n")
    
    else:
        print("--> FOUND all taxonomic files in", dirOut, "!\n")
    
    
    # BUILDING INDEX
    path_to_centriBuild = ("/home/sheldon/Applications/centrifuge-master21" +
                               "jan19/build/bin/centrifuge-build")
    db_name = "ncbi_16s"
    cmd_centriBuild = (path_to_centriBuild + " -p " + NB_THREADS + 
                       " --conversion-table " + dirOut + "seqid2taxid " + 
                       "--taxonomy-tree " + dirOut + 
                       "nodes.dmp --name-table " + dirOut + "names.dmp " + 
                       input_db_path + " " + dirOut + db_name)
    print("Building custom Centrifuge index...")
    # print(cmd_centriBuild);sys.exit()
    with open(dirOut + "centri-build.log", 'w') as centriBuild_log:
        sub.Popen(cmd_centriBuild.split(), stdout=centriBuild_log).communicate()
        
    print("Index building finished !")
    print("A log file can be found in", dirOut, '\n')
    