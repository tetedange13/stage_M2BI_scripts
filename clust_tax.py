#!/usr/bin/env python3

"""
Make the following steps from a given cluster of reads:

1) Pick randomly 1 read from the cluster (the 1st one acutally)
2) BLASTn this read against the nt-NCBI database [accesssion: 12 feb 2019]
3) Get the taxid associated with the top BLAST hit
4) Get the taxonomic rank associated with the taxid
5) Generate a report (csv file)

Usage:
  clust_tax.py (-i <inputFqFile>) (-c <clustFile>) [-m <minMemb>] [-a <altHits>] [-t <threads>] [-l <taxoCut>]
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inputFqFile=input_fq  input fastq file (given for clustering)
  -c --clustFile=clust       input txt file detailling the clusters
  -m --minMemb=min_memb      minimum number of members within a cluster [default: 250]
  -a --altHits=alt_hits      number of alternative BLAST hits to look for [default: 0]
  -l --taxoCut=taxo_cutoff   cutoff for the taxonomic level [default: none]
  -t --threads=nb_threads    number of threads to use [default: 10]
  
Remark: The  cutoff for the taxonomic level becomes mandatory when more than 1
        alternative BLAST hit is specified
"""

import sys, os, re
import Bio
import pandas as pd
import time as t
from functools import partial
import multiprocessing as mp
from Bio import Entrez, SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from docopt import docopt
import src.check_args as check
import src.parallelized as pll


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
    
    # Regex for taxid:
    regex_taxid = re.compile("(:?\/db_xref=\"taxon:)([0-9]+)") 
    # Regex for the ScientificName of the organism: 
    regex_organism = re.compile("(:?ORGANISM\s+)(.*)")
    # Not working on multiline taxonomy:
    regex_taxo = re.compile("(:?ORGANISM.*\n\s+)(.*)") 
    
    match_taxid = regex_taxid.search(fetch_res)
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
    

def dict_from_BLASTres(blast_res, idx_to_take):
    """
    Extract the different needed informations from a BLAST result and put them
    into a dictionnary, with exact same keys as the later dataFrame
    """
    dict_res = { "readID":blast_res.query.split()[0] }
    
    if len(blast_res.descriptions) == 0: # NO BLAST hit
        dict_res["problems"] = "NO_BLAST_HIT"
        return dict_res
        
    descr_blast = blast_res.descriptions[idx_to_take]
    if idx_to_take == 0:
        print("BEST:", descr_blast)
        print("bit-score (best HSP) =", 
              blast_res.alignments[idx_to_take].hsps[0].bits)
    
    taxid_hit, sp_name_hit, taxo_hit = annot_from_title(descr_blast.title)   
    
    dict_res["remarks"] = str(descr_blast)
    dict_res["topHit"] = sp_name_hit
    dict_res["taxid"] = taxid_hit
    dict_res["taxo"] = taxo_hit
    
    return dict_res


def query_rank_remote(str_taxid):
    """
    Given a taxid fetch the taxonomic rank querying (remotely) the Taxonomy db
    """
    print("SALUTUU", str_taxid)
    fetch_handle = Entrez.efetch(db="Taxonomy", id=str_taxid)
    res_fetch = Entrez.read(fetch_handle)
    fetch_handle.close()

    return res_fetch[0]["Rank"]
    

def query_rank_local(str_taxid):
    """
    Given a list of taxids, search their associated taxonomic ranks locally, by
    using the NCBI dump files (names and nodes) stored
    
    ENCORE EN CHANTIER...
    """
    return taxfoo.get_taxid_rank(str_taxid)


def look_for_alt(tupl_blast_res, nb_alt, my_df):
    """
    Take a BLAST result that has been identified as a FP (assigned but not 
    within the Zymo) and query other possible hits (for the same read), looking
    for hits with exact same score (and  e-value)
    """
    clust_key, blast_res = tupl_blast_res
    descr_topHit = blast_res.descriptions[0]
    splitted_remark = str(descr_topHit).split()
    # Reference values (from the top hit):
    max_e_val, max_score = splitted_remark[-1], splitted_remark[-2]
    topHit = descr_topHit.title.split(',')[0].split('| ')[1]
    
    print("\nTOP HIT:", topHit, " | ", "E-VAL:", max_e_val, " | ", "SC:", 
          max_score)
    
    alt = False
    for j in range(1, nb_alt+1):
        if len(blast_res.descriptions) > j: # If enough hits
            descr_alt = blast_res.descriptions[j]
            alt_score, alt_e_val = str(descr_alt.score), str(descr_alt.e)
            print("Looking for alternative hit #" + str(j))
            print(descr_alt)
            
            if (alt_score == max_score and alt_e_val == max_e_val):
                alt = True
                to_df_alt = dict_from_BLASTres(blast_res, j)
                keys_alt = to_df_alt.keys()
                alt_idx = "alt_" + clust_key.split('_')[1] + "_" + str(j)
                
                if "problems" in keys_alt:
                    print("PROBLEM WITH ALT!", to_df["problems"])
                
                else:
                    print("FOUND:", to_df_alt["topHit"], " | ", "E-VAL:", 
                    alt_e_val, " | ", "SC:", alt_score)
                    for key_alt in keys_alt:
                        my_df.loc[alt_idx, key_alt] = to_df_alt[key_alt]
                
    
    if not alt:
        print("NO acceptable alternative hit found")
    
   
def handle_strain(sp_name, rank):
    """
    Deals with 'strain' issues (i.e. organism that have been detected beyond
    the species level taxonomicly = 'no rank')
    Cut the name of the organism to keep only the 'Genus species' form and
    attributes the taxonomic level 'strain' to it  
    """
    splitted_sp_name = sp_name.split()
    
    # Very likely a strain if name longer than 2:
    if rank in ('no rank', 'subspecies') and len(splitted_sp_name) > 2: 
        # We keep only the first 2 words and we change the rank value:
        new_sp_name, new_rank_sp = " ".join(splitted_sp_name[0:2]), "strain"
        
        return (new_sp_name, new_rank_sp)
    
    else:
        return (sp_name, rank)
    


# MAIN:
if __name__ == "__main__":
    # COMMON VARIABLES AND PATHES:
    taxo_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 
                   'genus', 'species', 'subspecies'] + ["strain"]   
    ARGS = docopt(__doc__, version='0.1')
    input_fq_path, input_fq_base = check.infile(ARGS["--inputFqFile"],
                                                ('fastq', 'fq'))[0:2]
    TO_CLUST_FILE = ARGS["--clustFile"]
    NB_THREADS = check.input_nb(ARGS["--threads"])
    NB_MIN_BY_CLUST = check.input_nb(ARGS["--minMemb"]) # Needed later
    NB_ALT = int(check.input_nb(ARGS["--altHits"]))
    taxonomic_level_cutoff = check.acceptable_str(ARGS["--taxoCut"], 
                                                  taxo_levels + ['none'])
    
    # CHECK ARGS FOR ALTERNATIVE HITS AND CUTOFF:
    if NB_ALT == 1 and taxonomic_level_cutoff == "none":
        print("ERROR! As you want alternative hits, you have to specify a "
              "cutoff for the taxonomic level\n")
        sys.exit(2)
    elif NB_ALT > 1:
        print("Sorry, for now this script does not work with more than 1 " +
              "alternative BLAST hit (it largely enough most of the time)\n")
        sys.exit(2)
    
    out_xml_file = input_fq_base + "_memb" + NB_MIN_BY_CLUST + ".xml" # Needed later
    path_apps = "/home/sheldon/Applications/"
    START_TIME = t.time()
    
    if not os.path.isfile(out_xml_file):
        # GENERATE 1 FASTA FILE PER CLUSTER (with CARNAC-LR's dedicated script):
        to_CARNAC_to_fasta = (path_apps + "CARNAC-LR_git93dd640/scripts/" +
                              "CARNAC_to_fasta.py")
        out_clust_fa = "tmp_clust_memb" + NB_MIN_BY_CLUST + "/"
        
        cmd_clust_to_fasta = ("python3 " + to_CARNAC_to_fasta + " " +
                              TO_CLUST_FILE + " " + input_fq_path + " " + 
                              out_clust_fa + " " + NB_MIN_BY_CLUST)

        print("Generating 1 fasta file by cluster (with minimum", 
              NB_MIN_BY_CLUST, "members) ...")
        os.system(cmd_clust_to_fasta)
        
        # GET 1 READ BY CLUSTER AND CONCATENATE ALL THEM (into the stream):
        BLAST_TIME = t.time()
        cmd_pick_one = "head -qn 2 " + out_clust_fa + "cluster*.fasta"
        
        # QUERY THIS FASTA FILE TO THE NCBI (nt) DATABASE:
        to_blast_exe = path_apps + "blast-2.8.1+-src/ReleaseMT/bin/blastn"
        to_NCBIdb = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
                     "metage_ONT_2019/nt_db/BLASTdb/nt")
        
        cmd_BLAST = str(NcbiblastxCommandline(cmd=to_blast_exe, db=to_NCBIdb, 
                                              num_threads=NB_THREADS,
                                              out=out_xml_file, outfmt=5,
                                              max_target_seqs=100,
                                              word_size=11))
        cmd_BLAST += " -task=blastn" # Pb with this specific parameter
        #print(cmd_BLAST); sys.exit()
        print("BLAST query against NCBI db for 1 read by cluster...")
        os.system(cmd_pick_one + " | " + str(cmd_BLAST))
        
        #os.system("rm -r " + out_clust_fa) # Cleaning
        print("BLAST query finished !")
        print("TOTAL TIME FOR THE BLAST QUERY:", t.time() - BLAST_TIME, '\n')
    
    
    # PARSING BLAST OUTPUT, TO EXTRACT TAXID:
    Entrez.email = "felix.deslacs@gmail.com" # Config Entrez
    report_filename = input_fq_base + "_memb" + NB_MIN_BY_CLUST + ".csv"
    pb_filename = input_fq_base + "_memb" + NB_MIN_BY_CLUST + ".log"
    
    if not os.path.isfile(report_filename):
        print("PARSING BLAST OUTPUT...")

        PARSING_TIME = t.time()
        df_hits = pd.DataFrame(data=None)
        res_handle = open(out_xml_file)
        blast_records = NCBIXML.parse(res_handle)
        
        if os.path.isfile(pb_filename):
            os.remove(pb_filename)
        
        # GET TAXIDs AND OTHER INFORMATION:    
        for idx, blast_res in enumerate(blast_records):
            print("\n\nCLUSTER", idx, " | ", "QUERY =", 
                  blast_res.query.split()[0])
            
            to_df = dict_from_BLASTres(blast_res, 0) # First hit first
            keys_to_df = to_df.keys()
            
            if "problems" in keys_to_df: # Something went wrong somewhere
                print("PROBLEM! -->", to_df["problems"])
                with open(pb_filename, 'a') as pb_log_file:
                    pb_log_file.write("clust_" + str(idx) + " | " + 
                                      "QUERY: " + to_df["readID"] + " | " + 
                                      to_df["problems"] + '\n')
            
            else:
                # We can store informations about the top hit:
                keys_to_df = to_df.keys()
                idx_df = "clust_" + str(idx)
                
                for key_to_df in keys_to_df:
                    df_hits.loc[idx_df, key_to_df] = to_df[key_to_df]
        
        res_handle.close()
        
        # GET TAXONOMIC RANKS ASSOCIATED WITH EACH FOUND TAXID:
        query_rk = query_rank_remote # Pointer to function
        local_taxo_search = True
        
        if local_taxo_search:
            print("Requesting taxonomic ranks locally...")
            import src.ncbi_taxdump_utils as taxo_utils

            nodes_path = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019" + 
                          "/nt_db/taxo_18feb19/nodes.dmp")
            names_path = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019" + 
                          "/nt_db/taxo_18feb19/names.dmp")
            global taxfoo
            taxfoo = taxo_utils.NCBI_TaxonomyFoo()
            taxfoo.load_nodes_dmp(nodes_path)
            taxfoo.load_names_dmp(names_path)
            query_rk = query_rank_local # Pointer to function changed   
        
            # Parallelized local taxo rank search:
            partial_rk_search = partial(pll.rk_search, 
                                        query_rank_func=query_rk) 
            pool_parsing = mp.Pool(int(NB_THREADS))
            list_ranks = pool_parsing.map(partial_rk_search, 
                                          enumerate(df_hits["taxid"]))
            pool_parsing.close()
            
        else:
            print("\nRequesting taxonomic ranks remotely...")
            # Serial remote taxo rank search:    
            list_ranks = []
            for tupl_taxids in enumerate(df_hits["taxid"]):
                list_ranks.append(pll.rk_search(tupl_taxids, query_rk))
            
        # Add it to the dataFrame:
        for tupl_rak in list_ranks:
            current_df_idx, current_rk = tupl_rak
            df_hits.loc[current_df_idx, "rank"] = current_rk
        
        # Results saved to a report file:
        df_hits.to_csv(report_filename, sep='\t')
        print("\nTOTAL TIME FOR PARSING THE BLAST RESULTS:", 
              t.time() - PARSING_TIME, '\n')
        
        
    else:
        print("FOUND REPORT FILE!\n")
        df_hits = pd.read_csv(report_filename, sep='\t', index_col=0)  
    
    
    # LOOKING FOR ALTERNATIVE BLAST HITS:
    if NB_ALT > 0:# and len(dict_problems['FP']) > 0: # If there are some FP
        print("Looking for alternative BLAST hits (with EXACT same " +
              "score) for FP...") 

        res_handle = open(out_xml_file)
        blast_records = NCBIXML.parse(res_handle)
        func_parallel = partial(pll.FP_search, my_df=df_hits,
                                taxo_cutoff=taxonomic_level_cutoff,
                                taxo_levels=taxo_levels)
                                
        pool_searchFP = mp.Pool(int(NB_THREADS))
        list_FP_tmp = pool_searchFP.map(func_parallel, enumerate(blast_records))
        pool_searchFP.close()
        res_handle.close() 
    
        list_FP = [tupl[1:3] for tupl in list_FP_tmp if tupl[0] == "FP"]
        for tuple_blast_res in list_FP:
            look_for_alt(tuple_blast_res, NB_ALT, df_hits)

        print("Search for alternative hits for FP finished !\n")     
        # We have to re-write the new report file, with problems solved:u
        df_hits.to_csv(report_filename, sep='\t')
    
    
    # DEAL WITH THE DIFFERENT PROBLEMS ENCOUNTERED:
    print("Solving the encountered problems...")
    if os.path.isfile(pb_filename): # Check if there are problems to solve
        with open(pb_filename, 'r') as pb_log_file:
            problems_list = pb_log_file.read().splitlines()
        
        dict_problems = {'NO_BLAST_HIT':[]}
        for problem in problems_list:
            splitted_problm = problem.split(" | ")
            problm_name = splitted_problm[-1]
            clust_name = splitted_problm[0]

            if problm_name == "NO_BLAST_HIT": # Plus emmerdant a gerer ca..
                dict_problems['NO_BLAST_HIT'].append(clust_name)
            
            else:
                print("UNKNOWN PROBLEM")
                sys.exit()
           
        print("Problems solved !")    
        # We have to re-write the new report file, with problems solved:
        df_hits.to_csv(report_filename, sep='\t')
    
    else:
        print("NO problems to solve")
    
    
    print("\n\nTOTAL RUNTIME:", t.time() - START_TIME)  
      
