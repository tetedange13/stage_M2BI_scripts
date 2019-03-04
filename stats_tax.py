#!/usr/bin/env python3


"""
Mapping statistics computation

Usage:
  main.py (-i <inFile>)
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inFile=input_file     input file
  
Remark: Big supposition is that the SAM file does not contain secondary alignments
"""


import sys, os
import os.path as osp
import pandas as pd
import time as t
import subprocess as sub
from docopt import docopt
from src.shared_fun import handle_strain, in_zymo
import src.check_args as check


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
    to_infile, _, _, _, ext_infile = check.infile(ARGS["--inFile"], 
                                                  ["csv", "tsv", "sam"])
    taxo_cutoff = "genus"

    # Common variables:
    to_apps = "/home/sheldon/Applications/"
    taxo_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 
                   'genus', 'species', 'subspecies'] + ["strain"] 

    
    SAM_MODE = False
    if ext_infile == ".sam":
      SAM_MODE = True

    if SAM_MODE:
      print("SAM MODE !\n")
      to_samtools = to_apps + "samtools_0.1.19/samtools"

      cmd_get_mapped = (to_samtools + " view -S -F 4 " + to_infile +
                        " | awk '{print $1,$2,$3,$5,$9}'")
      # print(cmd_get_mapped);sys.exit()
      # Run process and get
      to_report_csv = "salut.csv"
      with open(to_report_csv, 'w') as csv_to_write:
        sub.Popen(cmd_get_mapped.split(), stdout=csv_to_write).communicate()
      sys.exit()
      START_CSV_PARSING = t.time()
      with open(to_report_csv, 'r') as my_csv:
          dict_csv = {}

          for line in my_csv:
            splitted_line = line.rstrip('\n').split('\t')
            readID, SAM_flag, target, MAPQ = splitted_line
            #dict_csv[readID] = {"targets":[target], "MAPQ":[MAPQ]}
            dict_csv[readID] = {"targets":target, "MAPQ":MAPQ}
            
            if 2048 in decompose_SAMflag(SAM_flag): #or 256 in decompose_SAMflag(SAM_flag): # Suppose that NO secondary
              # We supposed that primary is already in the keys
              continue # We skîp supplementary alignments
              #dict_csv[readID]["targets"].append(target)
              #dict_csv[readID]["MAPQ"].append(MAPQ)

              #print(dict_csv) ; sys.exit()

      #df_hits = pd.read_csv(to_report_csv, sep='\t', index_col=0, header=None)
      #print(dict_csv) 
      print("PARSING TIME:", str(t.time() - START_CSV_PARSING))

      for readID in dict_csv:
        pass
        # if len(dict_csv[readID]["targets"]) # Chimeric alignment

      sys.exit()

    else:
      df_hits = pd.read_csv(to_report_csv, sep='\t', index_col=0, header=0)

    # CALCULATION OF STATISTICS:
    cutoff_e_val = 0.00001 # 10^-5
    dict_stats = {'TN':0, 'FN':0, 'TP':0, 'FP':0}
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
                          taxo_cutoff, taxo_levels)
            

            if res == "FP":
                # Look for alternative hits:
                print("Not within the Zymo --> Look for alternative hit!")

                alt_index = "alt_" + clust.split('_')[1] + "_1"
                if alt_index in rownames_df:
                    res_alt = in_zymo(df_hits.loc[alt_index, "topHit"], 
                                      df_hits.loc[alt_index, "rank"], 
                                      taxo_cutoff, taxo_levels)
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


