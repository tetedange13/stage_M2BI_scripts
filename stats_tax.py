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
import sklearn.metrics as skm


def dict_stats_to_vectors(dict_res):
  """
  Convert a 
  """
  vec_pred, vec_true = [], []

  for res in dict_res:
    if res == 'TP':
      vec_true.extend([True]*dict_res[res]) # Positive in reality
      vec_pred.extend([True]*dict_res[res]) # Well predicted positive
    elif res == 'TN':
      vec_true.extend([False]*dict_res[res]) # Negative in reality
      vec_pred.extend([False]*dict_res[res]) # Well predicted negative
    elif res == 'FN':
      vec_true.extend([True]*dict_res[res]) # Positive in reality
      vec_pred.extend([False]*dict_res[res]) # But predicted negative
    else: # if res == 'FP':
      vec_true.extend([False]*dict_res[res]) # Negative in reality
      vec_pred.extend([True]*dict_res[res]) # But predicted positive

  return (vec_true, vec_pred)


def calc_recall(dict_res):
  # Sensitivity, hit rate, recall, or true positive rate
  # TPR = TP/(TP+FN)
  # Specificity or true negative rate
  # TNR = TN/(TN+FP) 
  # Precision or positive predictive value
  # PPV = TP/(TP+FP)
  # Negative predictive value
  # NPV = TN/(TN+FN)
  # Fall out or false positive rate
  # FPR = FP/(FP+TN)
  # False negative rate
  # FNR = FN/(TP+FN)
  # False discovery rate
  # FDR = FP/(TP+FP)
  # Overall accuracy
  # ACC = (TP+TN)/(TP+FP+FN+TN)

  # return dict_res['TP']/(dict_res['TP'] + dict_res['FP'])
  return dict_res['TP']/(dict_res['TP'] + dict_res['FN'])



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
    elif "toSilva" in infile_base:
      print("DB GUESSED: SILVA")
      to_seqid2taxid = to_dbs + "Centri_idxes/SILVA/seqid2taxid"
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
      # print(input_samfile.references);sys.exit()

      my_set = set()
      for alignment in input_samfile.fetch(until_eof=True):        
        if alignment.is_supplementary or alignment.is_secondary:
            print(alignment.query_name);sys.exit()
            continue # We skip secondary or supplementary alignments

        else:
          compt_alignments += 1
          if alignment.is_unmapped:
            dict_stats['FN'] += 1 # Count unmapped = 'FN'

          else:

            if alignment.mapping_quality >= mapq_cutoff: #and alignment.template_length > une_certaine_valeur
              ref_name = alignment.reference_name
              my_set.add(ref_name)
              #if re_name in ref_name:
              #  print(ref_name)
              # ref_name = alignment.reference_name.split()[0]
              current_taxid = dict_seqid2taxid[ref_name]
              assert(current_taxid) # If taxid not 'None'
              sp_name = taxfoo.get_taxid_name(current_taxid)
              current_rank = taxfoo.get_taxid_rank(current_taxid)
              if "SC" in alignment.reference_name:
                pass
                #print(in_zymo(sp_name, current_rank, taxo_cutoff))
                #print(alignment.__str__());sys.exit()
              
              assert(sp_name) # Not 'None'
              assert(current_rank) # Not 'None'
              res = in_zymo(sp_name, current_rank, taxo_cutoff)
              dict_stats[res] += 1
            else:
              dict_stats['FN'] += 1

      input_samfile.close()
      print("SAM PARSING TIME:", str(t.time() - START_SAM_PARSING))
      sum_stats = sum(dict_stats.values())
      print(sorted(dict_stats.items()))
      print("TOT ALIGNMENTS:", compt_alignments, " | ", "SUM STATS:", 
            sum_stats)
      assert(sum_stats == compt_alignments)

      y_true, y_pred = dict_stats_to_vectors(dict_stats)
      print("SENSITIVITY:", round(skm.recall_score(y_true, y_pred), 4))
      print("PRECISION:", round(skm.precision_score(y_true, y_pred), 4))
      print("F1-SCORE:", round(skm.f1_score(y_true, y_pred), 4))
      # print(skm.classification_report(y_true, y_pred))
      print("MATTHEWS:", round(skm.matthews_corrcoef(y_true, y_pred), 4))
      sys.exit()
      

    # else:
    #   df_hits = pd.read_csv(to_report_csv, sep='\t', index_col=0, header=0)

    # CALCULATION OF STATISTICS:
    cutoff_e_val = 0.00001 # 10^-5
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


