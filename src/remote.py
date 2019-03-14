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
    try:
        fetch_handle = Entrez.efetch(db="nucleotide", id=gi_str,
                                     rettype='gb', retmode="text")
        fetch_res = fetch_handle.read()
        fetch_handle.close()

    except HTTPError:
        print("ERROR! HTTP 400 WITH GI:", gi_str)
        sys.exit()
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
