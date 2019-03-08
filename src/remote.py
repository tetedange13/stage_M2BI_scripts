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


Entrez.email = "felix.deslacs@gmail.com" # Config Entrez


def query_taxid(term_to_search):
    """
    Given a term to search the taxonomy db, this function return the taxid
    (if one can be found)
    """
    try:
        search_handle = Entrez.esearch(db="Taxonomy", term=term_to_search)
        res_search = Entrez.read(search_handle)
        search_handle.close()
    except HTTPError:
        t.sleep(20) # Wait if 502 HTTP error from NCBI
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
        
        
def taxid_from_gb_entry(gi_str):
    """
    Fetch the GenBank entry using the GI of a sequence, read the entry and 
    parse it with regex, to get the taxid and the ScientificName of the 
    organism
    If the taxid cannot be found, the "Taxonomy" db of NCBI is (e-)searched
    using the ScientificName, to get the taxid 
    """  
    fetch_handle = Entrez.efetch(db="nucleotide", id=gi_str,
                                 rettype='gb', retmode="text")
    fetch_res = fetch_handle.read()
    fetch_handle.close()
    if not fetch_res:
        print("GI NOT FOUND:", gi_str)
        sys.exit()
    
    # Regex for taxid:
    regex_taxid = re.compile("(:?\/db_xref=\"taxon:)([0-9]+)")
    match_taxid = regex_taxid.search(fetch_res)
    
    # If we got the ScientficName, but not the taxid:
    if match_taxid:
        taxid = match_taxid.group(2)
    
    else:
        print("Taxid NOT contented within this GenBank entry !") 
        print("   --> Requesting it on the NCBI Taxonomy db")
        
        regex_organism = re.compile("(:?ORGANISM\s+)(.*)")
        match_organism = regex_organism.search(fetch_res)
        assert(match_organism)
        
        organism = match_organism.group(2)
        assert(organism != "metagenome")
        taxid = query_taxid(organism)
    
    # Case where the taxid CANNOT be queried with the given ScientificName  
    if taxid == "TAXO_NOT_FOUND" or taxid == "SEVERAL_TAXIDS":
        return "PROBLEM TAXID SEARCH"
        
    return taxid


def solve(to_problems_file):
    """
    For potential use in another part of the project, as an imported func
    """
    Entrez.email = "felix.deslacs@gmail.com" # Config Entrez
    acc2taxid_file = open("acc2taxid", 'a')

    with open(to_problems_file, 'r') as problems:
        for line in problems:
            t.sleep(1) # Make sure not too many requests to NCBI
            acc_nb = line.rstrip('\n')

            try:
                taxid = taxid_from_gb_entry(acc_nb)
            except HTTPError:
                t.sleep(20) # Wait if 502 HTTP error from NCBI
                taxid = taxid_from_gb_entry(acc_nb)

            acc2taxid_file.write(acc_nb + '\t' + taxid + '\n')
            #sys.exit()

    acc2taxid_file.close()


# MAIN:
if __name__ == "__main__":
    solve("problematic_acc_numbers.txt")
    
    # Je le garde, au cas où:
    # to_acc2taxid = "acc2taxid"
    # dict_acc2taxid = {}

    # with open(to_acc2taxid) as acc2taxid_file:
    #     for line in acc2taxid_file:
    #         acc_nb, taxid = line.rstrip('\n').split('\t')
    #         dict_acc2taxid[acc_nb] = taxid

    # list_gi_path = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
    #                 "metage_ONT_2019/rrn_8feb19/seqIDs_GIs.txt")
    # to_seqid2taxid = "seqid2taxid"

    # with open(to_seqid2taxid, 'w') as seqid2taxid_file, \
    #      open(list_gi_path, 'r') as list_gi_file:
    #      for line in list_gi_file:
    #         seqid, acc_nb = line.rstrip('\n').split('\t')
    #         seqid2taxid_file.write(seqid + '\t' + dict_acc2taxid[acc_nb] + 
    #                                '\n')
