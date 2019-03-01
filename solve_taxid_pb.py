#!/usr/bin/env python3

"""
Aims to solve problems linked with local taxids mapping (i.e. not found taxids),
by simply requeting directly the NCBI databases, using Biopython

The script has been adapted, to possibly handle a quite big number of taxids
to map

"""


import os, sys, re
import time as t
from Bio import Entrez
from urllib.error import HTTPError


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


# MAIN:
if __name__ == "__main__":
    Entrez.email = "felix.deslacs@gmail.com" # Config Entrez
    acc2taxid_file = open("acc2taxid", 'a')

    with open("problematic_acc_numbers.txt", 'r') as problems:
        for line in problems:
            t.sleep(1) # Make sure not too many requests to NCBI
            acc_nb = line.rstrip('\n')

            try:
                taxid = annot_from_title(acc_nb)
            except HTTPError:
                t.sleep(20) # Wait if 502 HTTP error from NCBI
                taxid = annot_from_title(acc_nb)

            acc2taxid_file.write(acc_nb + '\t' + taxid + '\n')
            #sys.exit()

    acc2taxid_file.close()
