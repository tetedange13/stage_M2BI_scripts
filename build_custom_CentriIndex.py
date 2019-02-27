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

import sys, os, gzip
import os.path as osp
import subprocess as sub
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
            
            if gi not in accDict.keys():
                accDict[gi] = seq_id
            else:
                #print("ERROR with dict conversion\n")
                continue
                #sys.exit(2)
            
    return accDict

def chunker_list(big_list, nb_chunks):
    return (big_list[i::nb_chunks] for i in range(nb_chunks))
    
def cut_list_in_chunks(big_list, nb_chunks):
    """
    Yield successive n-sized chunks from l.
    
    From: https://stackoverflow.com/questions/312443/ \
    how-do-you-split-a-list-into-evenly-sized-chunks
    """
    len_big_list = 257333726
    for i in range(0, len_big_list, nb_chunks):
        yield big_list[i:i + nb_chunks] 
    
def test_pll(chunk_acc2gi, list_accessions):
    for line in chunk_acc2gi:
        _, acc_nb, taxid, _ = line.split('\t')
    
        if acc_nb in list_accessions:
            print(acc_nb, taxid)
            return (acc_nb, taxid)
    
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
        
        
def annot_from_title(title):
    """
    Query the 'nucleotide' db to get the 'annotations' section associated
    with the GenBank entry corresponding to the GI contained in the 'title'
    of a BLAST result 
    """
    type_id, gi_hit = title.split('|')[0:2]
    assert(type_id == 'gi') # If it is not the GI, we have a problem
   
    fetch_handle = Entrez.efetch(db="nucleotide", id=gi_hit,
                                 rettype='gb', retmode="text")
    fetch_res = fetch_handle.read()
    fetch_handle.close()
    
    regex_taxid = re.compile("(:?\/db_xref=\"taxon:)([0-9]+)") # Regex for taxid
    regex_organism = re.compile("(:?ORGANISM\s+)(.*)")
     # Needed to get the taxonomy, even on several lines
    #regex_taxo = re.compile(, flags=re.DOTALL)
    # Not working on multiline taxonomy:
    regex_taxo = re.compile("(:?ORGANISM.*\n\s+)(.*)") 
    
    match_taxid = regex_taxid.search(fetch_res)
    # Regex for the ScientificName of the organism: 
    match_organism = regex_organism.search(fetch_res)
    match_taxo = regex_taxo.search(fetch_res)
    
    # If we got the ScientficName, but not the taxid:
    if match_organism and not match_taxid:
        print("Taxid NOT contented within this GenBank entry !") 
        print("   --> Requesting it on the NCBI Taxonomy db")
        organism = match_organism.group(2)
        taxid = query_taxid(organism)
        taxonomy = match_taxo.group(2)
    
    elif match_taxid and match_organism:
        taxid = match_taxid.group(2)
        # We keep the whole ScientificName (not just "Genus species"):
        organism = match_organism.group(2) 
        taxonomy = match_taxo.group(2)
        
    else:
        print("PB REGEX")
        sys.exit()
    
    # Case where the taxid CANNOT be queried with the given ScientificName  
    if taxid == "TAXO_NOT_FOUND":
        print("The ScientificName doesn't correspond to any taxid")
        print("    --> Querying the taxid with only 'Genus species'")
        taxid = query_taxid(' '.join(organism.split()[0:2]))
        
    return (taxid, organism, taxonomy)
    
    
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
        id_type = "gi"
        
        if id_type == "sp_name":
            Entrez.email = "felix.deslacs@gmail.com" # Config Entrez
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
            to_gz_acc2taxid = (to_dbs + "nt_db/taxo_18feb19/" +
                              "nucl_gb.accession2taxid.gz")
            separ = ' ' # Change the seperator if infile has different format
            dict_acc = accList_to_dict(list_gi_path, separ)
            list_acc = sorted(dict_acc.keys())
            seqid2map_file = open(to_seqid2map, 'w')

            # Open gz file in "read-text" mode (option "rt"):    
            with gzip.open(to_gz_acc2taxid, 'rt') as acc2taxid_gz:
                acc2taxid_gz.readline() # Skip header
                salut = acc2taxid_gz.read().splitlines()
                #salut = [ (line.split('\t')[1] + '\t' + line.split('\t')[2]) for line in acc2taxid_gz]
                print("LIST LOADED")
                chunks_acc2taxid = chunker_list(salut, NB_THREADS)
                
                # Parallel version:
                print("TOTO")
                pll_getTaxid = partial(test_pll, list_accessions=list_acc)
                pool_getTaxid = mp.Pool(NB_THREADS)
                pool_getTaxid.map(pll_getTaxid, chunks_acc2taxid)
                pool_getTaxid.close()
                print("FINISHED")
                
                sys.exit()
            
                for toto in salut:
                    print(toto)
                    sys.exit()
            
                #print(acc2taxid_gz.readline())
                #print(acc2taxid_gz.readline())
                #sys.exit()
                #for line in acc2taxid_gz:
                    
                #    _, acc_nb, taxid, _ = line.split('\t')
                    #print(acc_nb) ; sys.exit()
                    
                #    if acc_nb in list_acc:
                #        seqid2map_file.write(dict_acc[acc_nb] + '\t' + taxid)
                #        print("FOUND!", acc_nb, dict_acc[acc_nb],taxid)
                #        sys.exit()
                #        break
            
            seqid2map_file.close()
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

