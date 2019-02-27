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
1) The "-a" parameter becomes mandatory if no file called 'seqid2taxid.map' 
(mapping sequences IDs with taxids) could be found in 'dirOut/'
2) This file should be formated like that: 'seqid\tGI'
"""

import sys, os, gzip, re
import os.path as osp
import subprocess as sub
import time as t
import multiprocessing as mp
from Bio import Entrez
from functools import partial
from Bio.Blast import NCBIXML
from docopt import docopt
import src.check_args as check


def accList_to_dict(to_accList, separator):
    """
    Convert a file into a dict {acc_number:sequence_ID}
    """
    accDict = {}
    
    with open(to_accList, 'r') as accList_file:
        for line in accList_file:
            seq_id, gi = line.rstrip('\n').split()
            
            # if gi not in accDict.keys(): # Handle duplicate GI..
            if True:
                accDict[seq_id] = gi
            else:
                #print("ERROR with dict conversion\n")
                continue
                #sys.exit(2)
            
    return accDict

    
def query_taxid(term_to_search):
    """
    Given a term to search the taxonomy db, this function return the taxid
    (if one can be found)
    """
    search_handle = Entrez.esearch(db="Taxonomy", term=term_to_search)
    res_search = Entrez.read(search_handle)
    search_handle.close()
    
    taxid = res_search["IdList"]
    nb_taxid = len(taxid)

    if nb_taxid == 0: # No results found for the "term" specified
        return "TAXO_NOT_FOUND"
    elif nb_taxid > 1: # More than 1 taxid ("term" not enough precise ?)
        return "SEVERAL_TAXIDS"    
    else:
        return taxid[0]
        
        
def annot_from_title(gi_str):
    """
    Query the 'nucleotide' db to get the 'annotations' section associated
    with the GenBank entry corresponding to the GI contained in the 'title'
    of a BLAST result 
    """  
    fetch_handle = Entrez.efetch(db="nucleotide", id=gi_str,
                                 rettype='gb', retmode="text")
    fetch_res = fetch_handle.read()
    fetch_handle.close()
    
    # Regex for taxid  :
    regex_taxid = re.compile("(:?\/db_xref=\"taxon:)([0-9]+)")  
    match_taxid = regex_taxid.search(fetch_res)
    
    # If we got the ScientficName, but not the taxid:
    if match_taxid:
        taxid = match_taxid.group(2)
    
    else:
        print("Taxid NOT contented within this GenBank entry !") 
        print("   --> Requesting it on the NCBI Taxonomy db")
        taxid = query_taxid(organism)
    
    # Case where the taxid CANNOT be queried with the given ScientificName  
    if taxid == "TAXO_NOT_FOUND" or taxid == "SEVERAL_TAXIDS":
        return "PROBLEM TAXID SEARCH"
        
    return taxid
    
    
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
                                 ('fa', 'fasta', 'fq', 'fastq'))[0]
    dirOut = ARGS["--dirOut"]                            
    list_gi_path = ARGS["--accList"]
    NB_THREADS = 15
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    #dirIn_db = "/projets/metage_ONT_2019/databases/Zymo_genomes-ZR160406/"
    
    #global seqid2taxid_file
    #seqid2taxid_file = dirOut + "seqid2taxid.map"
    to_seqid2map = osp.join(dirOut + "seqid2taxid.map")
    
    if not osp.isfile(to_seqid2map) and list_gi_path == "none":
        print("ERROR! As 'dirOut/seqid2taxid.map' could not be found, you " +
              "need to specify a file containing the gi of the species " + 
              "contained in your database")
        sys.exit(2)
    
    
    #if list_gi_path != "none": # GETTING TAXIDS FROM NAME OR GIs:
    else:
        Entrez.email = "felix.deslacs@gmail.com" # Config Entrez
        id_type = "gi"
        
        if id_type == "sp_name":
            with open(dirIn_db + "Zymo_spList.txt", 'r') as spList_file:
                list_sp = spList_file.read().splitlines()

            list_taxids = []
            for sp_fasta in list_sp:
                sp_noFasta = osp.splitext(sp_fasta)[0]
                print("Getting taxid of", sp_noFasta)
                
                search_handle = Entrez.esearch(db="taxonomy", term=sp_noFasta)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                taxid = search_results['IdList'][0]
                list_taxids.append(taxid)
                
                with open(dirIn_db + sp_fasta, 'r') as input_fa_file:
                    map_headers_to_taxid(dirOut, input_fa_file, taxid)

        elif id_type == "gi":
            separ = ' ' # Change the seperator if infile has different format
            #dict_acc = accList_to_dict(list_gi_path, separ)
            #list_acc = sorted(dict_acc.keys())
            #set_acc = dict_acc.values()
            TAXID_START_TIME = t.time
            dict_memory = {} # To remember already queried accession numbers
            #seqid2map_file = open(to_seqid2map, 'w')
            
            if osp.isfile(to_seqid2map):
                os.remove(to_seqid2map)
            
            with open(list_gi_path, 'r') as list_gi_file:
                #compt = 0
                for line in list_gi_file:
                    #compt += 1
                    seq_id, acc_nb = line.rstrip('\n').split(separ)
                    
                    if acc_nb not in dict_memory.keys():
                        print("QUERYING", acc_nb, "FOR TAXID")
                        taxid = annot_from_title(acc_nb)
                        dict_memory[acc_nb] = taxid
                    else:
                        print("ACCESSION NUMBER ALREADY QUERIED")
                        taxid = dict_memory[acc_nb]
                    
                    # In any case, we write the file
                    with open(to_seqid2map, 'a') as seqid2map_file:
                        seqid2map_file.write(seq_id + '\t' + taxid + '\n')
                    
                    #if compt > 10:
                    #    sys.exit()
            
            
            #seqid2map_file.close()
            print("TOTAL TAXID SEARCH RUNTIME:", t.time() - TAXID_START_TIME)
            sys.exit()
            
                
            
    
    
    set_taxids = set()
    with open(to_seqid2map, 'r') as taxids_file:
        for line in taxids_file:
            set_taxids.add(line.split()[1])
    list_taxids = list(set_taxids)

    
    # GETTING DUMP FILES:
    if not (osp.isfile(dirOut+"nodes.dmp") and osp.isfile(dirOut+"names.dmp")):
        print("\nGenerating taxonomic tree and names...")
        os.system("centrifuge-download -o " + dirOut + " -t " + 
                  ','.join(list_taxids) + " taxonomy")
        print("Done !\n")
    
    else:
        print("All taxonomic files have been found !\n")
    
    
    # BUILDING INDEX
    path_to_centriBuild = ("/home/sheldon/Applications/centrifuge-master21" +
                               "jan19/build/bin/centrifuge-build")
    cmd_centriBuild = (path_to_centriBuild + " -p 7 --conversion-table " + 
                       dirOut + "seqid2taxid.map " + "--taxonomy-tree " + 
                       dirOut + "nodes.dmp --name-table " + dirOut + 
                       "names.dmp " + dirIn_db + "Zymo_genomes-ZR160406.fa " + 
                       dirOut + "Zymo")
    print("Building centrifuge index...")
    with open(dirOut + "centri-build.log", 'w') as centriBuild_log:
        sub.Popen(cmd_centriBuild.split(), stdout=centriBuild_log).communicate()
        
    print("Index building finished !")
    print("A log file can be found in", dirOut, '\n')
    


# TRASHES:
#search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
#data = Entrez.read(search)
#search.close()
#lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['phylum']}
#for d in data[0]['LineageEx']:
#    print(d)
#print(lineage)

#print(search.read())

