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
        print("   --> TAXO_NOT_FOUND FOR:", term_to_search)
        return "TAXO_NOT_FOUND"
    elif nb_taxid > 1: # More than 1 taxid ("term" not enough precise ?)
        print("   --> SEVERAL_TAXIDS FOR:", term_to_search)
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
    dict_specials = {"CP006690":"1051650",
                     "CDRP01002148":"562",
                     "FOMU01000028":"321267",
                     "CP001098":"373903",
                     "LIGM01000002":"337330",
                     "MCIE01000002":"46170",
                     "JTBM01000001":"28901",
                     "LKUO01004057":"3659",
                     "KH326984":"32644", # TAXID 32644 = NCBI unidentified
                     "KH326984":"32644",
                     "FAOO01000006":"1643428",
                     "DQ219316":"6738", # Unclassified Anomura
                     "LQ104849":"715199", # Streptomyces sp. DSM 24056 not found
                     "LQ104848":"631015", # Streptosporangium sp. DSM 24060 not found
                     "AF052914":"859426", # ScientificName = 'Telestula sp.'
                     "D85647":"63190", # ScientificName = '[Rhizoctonia] zeae
                     "JQ977511":"13687", # Sphingomonas sp. Dca12 => Sphingomonas
                     "AY165769":"71484"} # Flabellina sp. VV-2003 => Flabellina
                     # "":"",
                     # "GU138662":"2052682", # Pythium ultimum => Globisporangium ultimum
                     # "AY598647":"1448053" # Pythium intermedium => Globisporangium intermedium
                     # "AY598657":"2052683", # Pythium ultimum var. ultimum => Globisporangium ultimum var. ultimum
                     # "MUIH01000001":"2548426", # Gammaproteobacteria included
                     # "MUHZ01000001":"2548426"} # Gammaproteobacteria included

    if gi_str in dict_specials.keys():
        print("SPECIAL:", gi_str)
        # We return also the gi_str, to avoid unpack errors
        return (dict_specials[gi_str], gi_str)

    print("QUERYING:", gi_str)
    try:
        fetch_handle = Entrez.efetch(db="nucleotide", id=gi_str,
                                     rettype='gb', retmode="text")
        fetch_res = fetch_handle.read()
        fetch_handle.close()
    except HTTPError:
        print("ERROR! HTTP 400 WITH GI:", gi_str)
        return "HTTP400"

    if not fetch_res:
        print("GI NOT FOUND:", gi_str)
        return "GI_NOT_FOUND"
    
    # Regex for taxid and SCientificName of the organism:
    regex_taxid = re.compile("(:?\/db_xref=\"taxon:)([0-9]+)")
    match_taxid = regex_taxid.search(fetch_res)
    regex_organism = re.compile("(:?ORGANISM\s+)(.*)")
    match_organism = regex_organism.search(fetch_res)
    assert(match_organism)
    organism = match_organism.group(2)
    assert(organism != "metagenome")
    
    # If we got the ScientficName, but not the taxid:
    if match_taxid:
        taxid = match_taxid.group(2)
    else:
        print("Taxid NOT contented within this GenBank entry !") 
        print("   --> Requesting it on the NCBI Taxonomy db")
        taxid = query_taxid(organism)
        if taxid in ("TAXO_NOT_FOUND", "SEVERAL_TAXIDS"):
            print("   --> FAILED !")
        else:
            print("   --> SUCCESS !")

    return (taxid, organism)
