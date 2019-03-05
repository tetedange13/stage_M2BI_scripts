#!/usr/bin/env python3


"""
Mapping statistics computation

Usage:
  main.py (-i <inFile>) (-l <taxoCut>) [-c <numCut>]
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inFile=input_file     input file
  -l --taxoCut=taxo_cutoff   cutoff for the taxonomic level [default: none]
  -c --numCut=num_cutoff     cutoff on MAPQ (if SAM mode) or e-value (else) [default: none]
"""


import sys, os
import os.path as osp
import pandas as pd
import time as t
import subprocess as sub
import pysam as pys
from docopt import docopt
from src.shared_fun import handle_strain, in_zymo
import src.check_args as check
import src.shared_fun as shared


def to_full_name(header):
    converter = {"BS":"Bacillus subtilis", "SA":"Staphylococcus aureus",
                 "EF":"Enterococcus faecalis", "LB":"Lactobacillus fermentum", 
                 "LM":"Listeria monocytogenes",  "PA":"Pseudomonas aeruginosa", 
                 "SE":"Salmonella enterica", "EC":"Escherichia coli", 
                 "CN":"Cryptococcus neoformans", 
                 "SC":"Saccharomyces cerevisiae"}
    return converter[header[0:2]]
       

def decompose_SAMflag(str_SAMflag):
    """
    Inspired from: http://broadinstitute.github.io/picard/explain-flags.html
    """
    list_2powers = [2**idx for idx in range(12)]
    
    return [power for power in list_2powers if int(str_SAMflag) & power]
    
   

# MAIN:
if __name__ == "__main__":
    # ARGUMENTS:
    ARGS = docopt(__doc__, version='0.1')
    #to_report_csv = "cDNA_run1_sampl50k_memb5.csv"
    # infile_path = "test_16S_primLEN_mapped.csv"
    to_infile, infile_base, _, _, ext_infile = check.infile(ARGS["--inFile"], 
                                                         ["csv", "tsv", "sam"])
    num_cutoff = ARGS["--numCut"]
    taxo_cutoff = ARGS["--taxoCut"]

    # Common variables:
    to_apps = "/home/sheldon/Applications/"
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    # taxo_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 
    #                'genus', 'species', 'subspecies'] + ["strain"] 
    dict_stats = {'TN':0, 'FN':0, 'TP':0, 'FP':0}

    # Guess the database used:
    if "toZymo" in infile_base:
      print("DB GUESSED: Zymo")
      to_seqid2taxid = to_dbs + "Centri_idxes/Zymo/seqid2taxid"
    elif "toRrn" in infile_base:
      print("DB GUESSED: rrn")
      to_seqid2taxid = to_dbs + "Centri_idxes/rrn/seqid2taxid"
    else:
      print("Unkown database !\n")
      sys.exit(2)
    
    # Guess the "mode":
    SAM_MODE = False
    if ext_infile == ".sam":
      SAM_MODE = True

    if SAM_MODE:
      print("SAM MODE !\n")

      # Set MAPQ cutoff:
      if num_cutoff == 'none':
        mapq_cutoff = 0
      else:
        mapq_cutoff = int(check.input_nb(num_cutoff, "-c cutoff (MAPQ here)"))

      # To dump files:
      nodes_path = to_dbs + "nt_db/taxo_18feb19/nodes.dmp"
      names_path = to_dbs + "nt_db/taxo_18feb19/names.dmp"
      
      print("Loading taxonomic Python module...")
      import src.ncbi_taxdump_utils as taxo_utils
      taxfoo = taxo_utils.NCBI_TaxonomyFoo()
      taxfoo.load_nodes_dmp(nodes_path)
      taxfoo.load_names_dmp(names_path)
      print("Taxonomic Python module loaded !\n")


      # To make correspond operon number and taxid:
      dict_seqid2taxid = {}
      with open(to_seqid2taxid, 'r') as seqid2taxid_file:
        for line in seqid2taxid_file:
          splitted_line = line.rstrip('\n').split('\t')
          dict_seqid2taxid[splitted_line[0]] = int(splitted_line[1])

      # Start SAM file parsing:
      print("Extracting information from SAM file...")
      START_SAM_PARSING = t.time()
      compt_alignments = 0
      input_samfile = pys.AlignmentFile(to_infile, "r")


      for alignment in input_samfile.fetch(until_eof=True):
        if alignment.is_unmapped:
          dict_stats['FN'] += 1 # Count unmapped = 'FN'
        else:
          if alignment.is_supplementary or alignment.is_secondary:
            continue # We skip secondary or supplementary alignments

          else:
            compt_alignments += 1
            if alignment.mapping_quality >= mapq_cutoff:
              current_taxid = dict_seqid2taxid[alignment.reference_name]
              assert(current_taxid) # If taxid not 'None'
              sp_name = taxfoo.get_taxid_name(current_taxid)
              current_rank = taxfoo.get_taxid_rank(current_taxid)
              assert(sp_name) # Not 'None'
              assert(current_rank) # Not 'None'
              res = in_zymo(sp_name, current_rank, taxo_cutoff)
              dict_stats[res] += 1
            else:
              dict_stats['FN'] += 1

      input_samfile.close()
      print("SAM PARSING TIME:", str(t.time() - START_SAM_PARSING))
      sum_stats = sum(dict_stats.values())
      print(dict_stats)
      print("TOT ALIGNMENTS:", compt_alignments, " | ", "SUM STATS:", 
            sum_stats)
      assert(sum_stats == compt_alignments)
      sys.exit()

      to_samtools = to_apps + "samtools_0.1.19/samtools"

      # cmd_get_mapped = (to_samtools + " view -S -F 4 " + to_infile +
      #                   " | awk '{print $1,$2,$3,$5,$9}'")

      to_report_csv = "tmp.csv"
      if not osp.isfile(to_report_csv):
        cmd_get_mapped = to_samtools + " view -S -F 4 " + to_infile
        cmd_count_unmapped = (to_samtools + " view -S -f 4 " + to_infile + 
                              " | wc -l")

        # Run process and get
        
        print("Writting new SAM file with mapped reads only and count",
              "unmapped reads...")
        START_WRITTING = t.time()
        with open(to_report_csv, 'w') as csv_to_write:
          sub.Popen(cmd_get_mapped.split(), stdout=csv_to_write).communicate()
        
        nb_unmapped = sub.Popen(cmd_count_unmapped.split(), 
                                stdout=sub.PIPE).communicate()[0]
        print()
        print("WRITTING AND COUNTING TIME:", str(t.time() - START_WRITTING))

      # EXTRACT INFO FROM SAM:
      print("Extracting information from SAM file...")
      START_CSV_PARSING = t.time()
      with open(to_report_csv, 'r') as my_csv:
          dict_csv = {}

          for line in my_csv:
            splitted_line = line.rstrip('\n').split('\t')[0:9]
            readID, SAM_flag, target, _, MAPQ, _, _, _, tlen = splitted_line
            #dict_csv[readID] = {"targets":[target], "MAPQ":[MAPQ]}
            if readID not in dict_csv.keys():
              dict_csv[readID] = {"targets":target, "MAPQ":MAPQ, 
                                  "len_align":tlen}
            
            
            else:
              if 2048 in decompose_SAMflag(SAM_flag): #or 256 in decompose_SAMflag(SAM_flag): # Suppose that NO secondary
              # We suppose that primary is already in the keys
                continue # We skîp supplementary alignments
                #dict_csv[readID]["targets"].append(target)
                #dict_csv[readID]["MAPQ"].append(MAPQ)

      #df_hits = pd.read_csv(to_report_csv, sep='\t', index_col=0, header=None)
      #print(dict_csv)  
      print("PARSING TIME:", str(t.time() - START_CSV_PARSING))
      os.remove(to_report_csv)
      sys.exit()

      # for readID in dict_csv:
        # pass
        # if len(dict_csv[readID]["targets"]) # Chimeric alignment

      

    else:
      df_hits = pd.read_csv(to_report_csv, sep='\t', index_col=0, header=0)

    # CALCULATION OF STATISTICS:
    cutoff_e_val = 0.00001 # 10^-5
    # dict_stats = {'TN':0, 'FN':0, 'TP':0, 'FP':0}
    dict_str = {'TN': "Assigned et un de nos Euk",
                'TP': "Assigned et une de nos bactos !",
                'FN': "Not assigned ! (cutoff plus bas dans l'arbre)",
                'FP': "Assigned, mais appartient pas à la Zymo"}
    
    rownames_df = df_hits.index
    rows_clust = [rowname for rowname in rownames_df if "clust" in rowname]

    for clust in rows_clust:
        splitted_remarks = df_hits.loc[clust, "remarks"].split()
        e_val  = splitted_remarks[-1]
        score = splitted_remarks[-2]
        
        print("\n", clust, " | ", "E-VAL:", e_val, " | ", "SCORE:", score)   
        print("NAME:", df_hits.loc[clust, "topHit"], " | ", 
              "RANK:", df_hits.loc[clust, "rank"])
        
        if float(e_val) < cutoff_e_val:
            res = in_zymo(df_hits.loc[clust, "topHit"], 
                          df_hits.loc[clust, "rank"], 
                          taxo_cutoff)
            

            if res == "FP":
                # Look for alternative hits:
                print("Not within the Zymo --> Look for alternative hit!")

                alt_index = "alt_" + clust.split('_')[1] + "_1"
                if alt_index in rownames_df:
                    res_alt = in_zymo(df_hits.loc[alt_index, "topHit"], 
                                      df_hits.loc[alt_index, "rank"], 
                                      taxo_cutoff)
                    if res_alt != 'FP':
                        print("FOUND:", df_hits.loc[alt_index, "topHit"], 
                              " | ", "RANK:", df_hits.loc[alt_index, "rank"]) 
                        remarks_alt = df_hits.loc[alt_index, "remarks"]
                        splitted_remarks_alt = remarks_alt.split()
                        print(alt_index, " | ", "E-VAL:", 
                              splitted_remarks_alt[-1], " | ", "SCORE:", 
                              splitted_remarks_alt[-2])        
                    else:
                        print("The alternative is still a FP")
                        
                    dict_stats[res_alt] += 1
                    print(dict_str[res_alt])
                        
                else:
                    print("NO alternative hit available")
                    dict_stats[res] += 1 
                    
            
            else:
                print(dict_str[res])
                dict_stats[res] += 1    
        
        
        else: # Not assigned because of the e_val cutoff
            print("Not assigned because of the e_val cutoff")
            dict_stats['FN'] += 1     
    
    print(dict_stats)


