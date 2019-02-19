#!/usr/bin/env python3

import sys, os
import subprocess as sub
from Bio import Entrez
    

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
    Entrez.email = "felix.deslacs@gmail.com" # Config Entrez
    #path_to_sh_script = "/projets/metage_ONT_2019/test.sh"
    dirIn_db = "/projets/metage_ONT_2019/databases/Zymo_genomes-ZR160406/"
    dirOut = "centri_idx-Zymo/"
    global seqid2taxid_file
    seqid2taxid_file = dirOut + "seqid2taxid.map"
    
    if not os.path.isdir(dirOut):
            os.mkdir(dirOut)
    
    
    # GETTING TAXIDS:
    if not os.path.isfile(seqid2taxid_file):
        with open(dirIn_db + "Zymo_spList.txt", 'r') as spList_file:
            list_sp = spList_file.read().splitlines()

        list_taxids = []
        for sp_fasta in list_sp:
            sp_noFasta = os.path.splitext(sp_fasta)[0]
            print("Getting taxid of", sp_noFasta)
            
            search_handle = Entrez.esearch(db="taxonomy", term=sp_noFasta)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            taxid = search_results['IdList'][0]
            list_taxids.append(taxid)
            
            with open(dirIn_db + sp_fasta, 'r') as input_fa_file:
                map_headers_to_taxid(dirOut, input_fa_file, taxid)

    else:
        set_taxids = set()
        with open(seqid2taxid_file, 'r') as taxids_file:
            for line in taxids_file:
                set_taxids.add(line.split()[1])
        list_taxids = list(set_taxids)

    
    # GETTING DUMP FILES:
    if not (os.path.isfile(dirOut+"nodes.dmp") and os.path.isfile(dirOut+"names.dmp")):
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

