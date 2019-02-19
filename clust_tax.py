#!/usr/bin/env python3

"""
Make the following steps from a given cluster of reads:

1) Pick randomly 1 read from the cluster (the 1st one acutally)
2) BLASTn this read against the nt-NCBI database [accesssion: 12 feb 2019]
3) Get the taxid associated with the top BLAST hit
4) Get the taxonomic rank associated with the taxid
5) Generate a report (csv file)

Usage:
  clust_tax.py (-i <inputFqFile>) (-c <clustFile>) [-m <minMemb>] [-a <altHits>] [-t <threads>]
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inputFqFile=input_fq  input fastq file (given for clustering)
  -c --clustFile=clust       input txt file detailling the clusters
  -m --minMemb=min_memb      minimum number of members within a cluster [default: 5]
  -a --altHits=alt_hits      number of alternative BLAST hits to look for [default: 1]
  -t --threads=nb_threads    number of threads to use [default: 10]
"""

import sys, os, re
import Bio
import pandas as pd
import time as t
import multiprocessing as mp
from Bio import Entrez, SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from docopt import docopt


def check_input_fq(input_fq_path):
    """
    Check if everything's fine with the fq file given as input
    """
    if not os.path.isfile(input_fq_path):
        print("ERROR! Wrong specified input fastq file")
        sys.exit(2)
    else:
        head_input_fq, tail_input_fq = os.path.split(input_fq_path)
        input_fq_base, input_fq_ext = os.path.splitext(tail_input_fq)
        
        if input_fq_ext[1:] not in ("fq", "fastq"):
            err_ext = ("ERROR! Current only fastq file (ext=fq or " +
                       "ext=fastq) are accepted")
            print(err_ext)
            sys.exit(2)
        
        else:
            return input_fq_path


def check_input_nb(input_nb):
    """
    Check if a integer given as argument is valid
    """
    try:
        int(input_nb)
        
    except ValueError:
        print("ERROR! You should give an integer here")
        sys.exit(2)
    
    return input_nb


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
        #organism = " ".join(match_taxo_org.group(2).split()[0:2])
        # We keep the whole name (not just "Genus species"):
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
    #else:
    #    print("Taxid successfully got")
        
    return (taxid, organism, taxonomy)


def query_taxonomy(str_taxid):
    """
    Given a taxid to fetch from the taxonomy db, this function returns the 
    the rank associated with the given taxid
    """
    fetch_handle = Entrez.efetch(db="Taxonomy", id=str_taxid)
    res_fetch = Entrez.read(fetch_handle)
    fetch_handle.close()
    lin_ex = res_fetch[0]["LineageEx"]

    return res_fetch[0]
    
    
def dl_gbk_file():
    """
    """
    # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&retmode=text&rettype=gb&id=835164530&email=felix.deslacs%40gmail.com
    # --> Est-ce que c'est vraiment une bonne idee ? Je vais encore me prendre
    # des stops par le site si je dl trop de trucs en meme temps...
    

def dict_from_BLASTres(blast_res, idx_to_take):
    """
    Extract the different needed informations from a BLAST result and put them
    into a dictionnary, with exact same keys as the later dataFrame
    """
    dict_res = { "readID":blast_res.query.split()[0] }
    
    if len(blast_res.descriptions) == 0: # NO BLAST hit
        #print("No BLAST hits --> ignored ?")
        dict_res["problems"] = "NO_BLAST_HIT"
        return dict_res
        
    descr_blast = blast_res.descriptions[idx_to_take]
    if idx_to_take == 0:
        print("BEST:", descr_blast)
        print("bit-score (best HSP) =", 
              blast_res.alignments[idx_to_take].hsps[0].bits)
    
    taxid_hit, sp_name_hit, taxo_hit = annot_from_title(descr_blast.title)
    
    #if taxid_hit == "TAXO_NOT_FOUND":
    #    dict_res["problems"] = taxid_hit
    #    return dict_res
    
    tax_record = query_taxonomy(taxid_hit)
    dict_res["remarks"] = str(descr_blast)
    dict_res["topHit"] = sp_name_hit
    dict_res["taxo"] = taxo_hit
    dict_res["rank"] = tax_record['Rank']

    return dict_res
    

def parse_BLAST(enum_iter):
    """
    Tiny function, mostly useful for parallelization actually
    """
    idx, blast_res = enum_iter
    
    readID = blast_res.query.split()[0]
    print("\n\nCLUSTER", idx, " | ", "QUERY =", readID)
    
    return ( idx, dict_from_BLASTres(blast_res, 0) )


def look_for_alt(blast_res, nb_alt, cutoff_e_val, my_df):
    """
    Take a BLAST result that has been identified as a FP (assigned but not 
    within the Zymo) and query other possible hits (for the same read), looking
    for hits with exact same score (and  e-value)
    """
    descr_topHit = blast_res.descriptions[0]
    splitted_remark = str(descr_topHit).split()
    # Reference values (from the top hit):
    max_e_val, max_score = splitted_remark[-1], splitted_remark[-2]
    topHit = descr_topHit.title.split(',')[0].split('| ')[1]
    
    print("\nTOP HIT:", topHit, " | ", "E-VAL:", max_e_val, " | ", "SC:", 
          max_score)
    
    if (float(max_e_val) < cutoff_e_val): 
        print("Alternatives hits (with EXACT SAME SCORE):")
        alt = False
        
        for j in range(1, nb_alt+1):
        #if False:
            descr_alt = blast_res.descriptions[j]
            alt_score, alt_e_val = str(descr_alt.score), str(descr_alt.e)
            
            if (alt_score == max_score and alt_e_val == max_e_val):
                alt = True
                to_df_alt = dict_from_BLASTres(blast_res, j)
                keys_alt = to_df_alt.keys()
                alt_idx = "alt_" + str(idx) + "_" + str(j)
                
                if "problems" in keys_alt:
                    print("PROBLEM WITH ALT!", to_df["problems"])
                    sys.exit()
                #    with open(pb_filename, 'a') as pb_log_file:
                #        pb_log_file.write(alt_idx + " | " + "QUERY =" +
                #                          to_df_alt + '\n')
                
                else:
                    print("FOUND:", to_df_alt["topHit"], " | ", "E-VAL:", 
                    alt_e_val, " | ", "SC:", alt_score)
                    for key_alt in keys_alt:
                        my_df.loc[alt_idx, key_alt] = to_df_alt[key_alt]
                
        
        if not alt:
            print("No alternative hits")
    
   
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


def in_zymo(sp_name, sp_rank, taxo_level_cutoff):
    """
    Given the rank of a BLAST hit and a sublist of taxonomic levels, deducted
    from the cutoff (for the taxonomic level),  returns a string in 
    ('TP', 'FP', 'TN', 'FN').
    The name of the species is also requiered, to deal with 'strain' cases
    """
    taxo_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 
                   'genus', 'species', 'subspecies'] + ["strain"]
    # Define the species contained within the Zymo mock community:
    list_prok = ['Listeria monocytogenes', 'Bacillus subtilis', 
                 'Staphylococcus aureus', 'Escherichia coli', 
                 'Lactobacillus fermentum', 'Enterococcus faecalis',
                 'Pseudomonas aeruginosa', 'Salmonella enterica']
    list_euk = ['Cryptococcus neoformans', 'Saccharomyces cerevisiae']
    sublist_taxo = taxo_levels[taxo_levels.index(taxo_level_cutoff): ]
    
    # We handle the issues associated with 'strains'
    new_name, rew_rank = handle_strain(sp_name, sp_rank)
    
    if rew_rank in sublist_taxo: # Assigned
        if new_name in list_prok:
            return 'TP' # True Positive
        elif new_name in list_euk:
            return 'TN' # True Negative
        else:
            return 'FP' # False Positive
    
    else:
        return 'FN'


# MAIN:
if __name__ == "__main__":
    # COMMON VARIABLES AND PATHES:
    path_apps = "/home/sheldon/Applications/"
    path_proj = "/projets/metage_ONT_2019/"
    out_xml_file = "../opuntia_" + NB_MIN_BY_CLUST + ".xml" # Needed later
    #NB_MIN_BY_CLUST = str(250) # 10 clusters
    #NB_MIN_BY_CLUST = str(25) # 20 clusters
    #NB_MIN_BY_CLUST = str(15) # 50 clusters
    
    ARGS = docopt(__doc__, version='0.1')
    input_fq_path = check_input_fq(ARGS["--inputFqFile"])
    CLUST_FILE = ARGS["--clustFile"]
    NB_THREADS = check_input_nb(ARGS["--threads"])
    NB_MIN_BY_CLUST = check_input_nb(ARGS["--minMemb"]) # Needed later
    NB_ALT = check_input_nb(ARGS["--altHits"])
    
    if not os.path.isfile(out_xml_file):
    #if True:         
        # GENERATE 1 FASTA FILE PER CLUSTER (with CARNAC-LR's dedicated script):
        to_CARNAC_to_fasta = (path_apps + "CARNAC-LR_git93dd640/scripts/" +
                              "CARNAC_to_fasta.py")
        to_CARNAC_txt = "../CARNAC/cDNA_run1_sampl50k_CARNAC.txt"
        to_fq = "../../cDNA_run1_sampl50k.fq"
        out_clust_fa = "tmp_clust/"
        
        cmd_clust_to_fasta = ("python3 " + to_CARNAC_to_fasta + " " +
                              to_CARNAC_txt + " " + to_fq + " " + out_clust_fa + 
                              " " + NB_MIN_BY_CLUST)
        #print(cmd_clust_to_fasta)
        print("Generating 1 fasta file by cluster (with minimum", 
              NB_MIN_BY_CLUST, "members) ...")
        os.system(cmd_clust_to_fasta)
        #sys.exit()
        
        # GET 1 READ BY CLUSTER AND CONCATENATE ALL THEM (into the stream):
        BLAST_TIME = t.time()
        #cmd_pick_one = ("head -qn 4 " + out_clust_fa + "cluster*.fasta | " +
        #                "sed -n '1~4s/^@/>/p;2~4p'")
        cmd_pick_one = "head -qn 2 " + out_clust_fa + "cluster*.fasta"
        #print(cmd_pick_one)
        
        # QUERY THIS FASTA FILE TO THE NCBI (nt) DATABASE:
        to_blast_exe = path_apps + "blast-2.8.1+-src/ReleaseMT/bin/blastn"
        to_NCBIdb = ("/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/" +
                     "metage_ONT_2019/nt_BLASTdb/nt")
        
        cmd_BLAST = NcbiblastxCommandline(cmd=to_blast_exe, db=to_NCBIdb, 
                                          num_threads=NB_THREADS,
                                          out=out_xml_file, outfmt=5)
        
        print("BLAST query against NCBI db for 1 read by cluster...")
        #print(cmd_BLAST) ; sys.exit()
        os.system(cmd_pick_one + " | " + str(cmd_BLAST))
        #stdout, stderr = toto() # Not needed cuz BLAST is not outputting anything
        
        #os.system("rm -r " + out_clust_fa) # Cleaning
        print("BLAST query finished !")
        print("TOTAL TIME FOR THE BLAST QUERY:", t.time() - BLAST_TIME, '\n')
    
    
    # PARSING BLAST OUTPUT, TO EXTRACT TAXID:
    Entrez.email = "felix.deslacs@gmail.com" # Config Entrez
    taxonomy_level_cutoff = "genus" 
    report_filename = "report_" + NB_MIN_BY_CLUST + ".csv"
    global pb_filename
    pb_filename = "problems_" + NB_MIN_BY_CLUST + ".log"
    
    if not os.path.isfile(report_filename):
    #if True:
        print("PARSING BLAST OUTPUT...")
        
        if NB_ALT > 0:
            list_FP = []
               
        PARSING_TIME = t.time()
        df_hits = pd.DataFrame(data=None)
        #NB_CLUSTERS = 0 # Count the total number of clusters
        res_handle = open(out_xml_file)
        #blast_records = SearchIO.parse(result_handle, 'blast-xml')
        blast_records = NCBIXML.parse(res_handle)
        
        if os.path.isfile(pb_filename):
            os.remove(pb_filename)
        
        # SERIAL VERSION:    
        for enum_iterable in enumerate(blast_records):
            #if enum_iterable[0] == 65:
            if True:
                current_blast_res = enum_iterable[1]
                idx, to_df = parse_BLAST(enum_iterable)
                keys_to_df = to_df.keys()
                
                if "problems" in keys_to_df: # Something went wrong somewhere
                    print("PROBLEM! -->", to_df["problems"])
                    with open(pb_filename, 'a') as pb_log_file:
                        pb_log_file.write("clust_" + str(idx) + " | " + 
                                          "QUERY: " + to_df["readID"] + " | " + 
                                          to_df["problems"] + '\n')
                
                else:
                    membership = in_zymo(to_df["topHit"], to_df["rank"], 
                                         taxonomy_level_cutoff)
                    if membership == 'FP' and NB_ALT > 0: 
                        list_FP.append(current_blast_res)
                        # We write it into the 'problems' file:
                        with open(pb_filename, 'a') as pb_log_file:
                            pb_log_file.write("clust_" + str(idx) + " | " + 
                                              "QUERY: " + to_df["readID"] + 
                                              " | " + membership + '\n')    
                    else:
                        # We can store informations about the top hit:
                        keys_to_df = to_df.keys()
                        idx_df = "clust_" + str(idx)
                        
                        for key_to_df in keys_to_df:
                            df_hits.loc[idx_df, key_to_df] = to_df[key_to_df]
        
        # PARALLELIZED VERSION:
        #pool_parsing = mp.Pool(int(NB_THREADS))
        #pool_parsing = mp.Pool(5)
        #pool_parsing.map(parse_BLAST, enumerate(blast_records), chunksize=8)
        #try:
        #    pool_parsing.map(parse_BLAST, enumerate(blast_records))
        #except mp.pool.MaybeEncodingError:
        #    print("Bonjour")
        
        #pool_parsing.close()
        #pool_parsing.join()

        res_handle.close()
        # Results saved to a report file:
        df_hits.to_csv(report_filename, sep='\t')
        print("\nTOTAL TIME FOR PARSING THE BLAST RESULTS:", 
              t.time() - PARSING_TIME, '\n')
        
        
    else:
        print("FOUND REPORT FILE!\n")
        df_hits = pd.read_csv(report_filename, sep='\t', index_col=0)  
    
    
    # DEAL WITH THE DIFFERENT PROBLEMS ENCOUNTERED:
    print("Solving the encountered problems...")
    with open(pb_filename, 'r') as pb_log_file:
        problems_list = pb_log_file.read().splitlines()
    
    dict_problems = {'FP':[], 'NO_BLAST_HIT':[]}
    for problem in problems_list:
        splitted_problm = problem.split(" | ")
        problm_name = splitted_problm[-1]
        clust_name = splitted_problm[0]
        
        if problm_name == 'FP':
            dict_problems['FP'].append(clust_name)
 
        elif problm_name == "NO_BLAST_HIT": # Plus emmerdant a gerer ca..
            dict_problems['NO_BLAST_HIT'].append(clust_name)
        
        else:
            print("UNKNOWN PROBLEM")
            sys.exit()
    
    if NB_ALT > 0 and len(dict_problems['FP']) > 0: # If there are some FP
        print("Looking for alternative BLAST hits for FP...") 
        if os.path.isfile(report_filename): # list_TP does not already exist
            res_handle = open(out_xml_file)
            blast_records = NCBIXML.parse(res_handle)
            list_FP = []
            for idx, blast_rec in enumerate(blast_records):
                if "clust_" + str(idx) in dict_problems['FP']:
                    list_FP.append(blast_rec)    
            res_handle.close() 
    
        for blast_res in list_FP:
            look_for_alt(blast_res, NB_ALT, 0.001, df_hits)
        
    print("Problems solved !")    
    # We have to re-write the new report file, with problems solved:
    df_hits.to_csv(report_filename, sep='\t')
    
