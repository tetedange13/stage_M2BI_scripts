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
import src.check_args as check
import sklearn.metrics as skm
import src.remote as remote
import matplotlib.pyplot as plt


def eval_align(ref_name, dict_conv_seqid2taxid, taxonomic_cutoff):
  """
  """
  current_taxid = dict_conv_seqid2taxid[ref_name]
  sp_name = taxfoo.get_taxid_name(current_taxid)
  
  if not sp_name:
    #print(ref_name, "HAS A PROBLEM !");sys.exit()
    to_problems_file = "problematic_ref_names.txt"
    if osp.isfile(to_problems_file):
      with open(to_problems_file, 'r') as problems:
        set_taxids = {line.rstrip('\n').split('\t')[1] for line in problems}
    else:
      set_taxids = set()

    if str(current_taxid) not in set_taxids:
      queried_taxid = remote.query_taxid(ref_name) # Try remotely
      if queried_taxid == "TAXO_NOT_FOUND": # Still not working
        with open(to_problems_file, 'a') as problems:
          problems.write(ref_name + '\t' + str(current_taxid) + '\n')
        print("PROBLEM with:", ref_name, "of associated taxid:", current_taxid)
        return None
      else:
        print("LE REMOTE A MARCHE !")
        sp_name = taxfoo.get_taxid_name(queried_taxid)

    else:
      return None

  assert(sp_name) # Not 'None'
  current_rank = taxfoo.get_taxid_rank(current_taxid)
  assert(current_rank) # Not 'None'

  return in_zymo(sp_name, current_rank, taxonomic_cutoff)


def handle_suppl(list_align_obj, dict_conv_seqid2taxid):
  """
  Take an alignment that has been found to be supplementary
  """
  set_taxid_target = set()
  for alignment in list_align_obj:
    taxid_target = dict_conv_seqid2taxid[alignment.reference_name]
    assert(taxid_target) # If taxid not 'None'
    set_taxid_target.add(taxid_target)
  
  # Different target name, but corresponding to the same taxid ?
  if len(set_taxid_target) == 1: 
    return alignment.reference_name # Any target name is OK
  return False

def dict_stats_to_vectors(dict_res):
  """
  Generate 2 vectors (for predicted and true values), based on the values of
  'TP', 'FP' etc contained in the dict_res given as arg
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


def calc_specificity(dict_res):
  # Negative predictive value: NPV = TN/(TN+FN)
  # FOR = 1 - NPV
  try:
    spe = dict_res['TN']/(dict_res['TN'] + dict_res['FP'])
    return spe
  except ZeroDivisionError:
    return 'NA'


def calc_FOR(dict_res):
  # Sensitivity, hit rate, recall, or true positive rate: TPR = TP/(TP+FN)
  # Specificity or true negative rate: TNR = TN/(TN+FP) 
  # Precision or positive predictive value: PPV = TP/(TP+FP)
  # Negative predictive value: NPV = TN/(TN+FN)
  # Fall out or false positive rate: FPR = FP/(FP+TN)
  # False negative rate: FNR = FN/(TP+FN)
  # False discovery rate: FDR = FP/(TP+FP)
  # Overall accuracy: ACC = (TP+TN)/(TP+FP+FN+TN)

  # return dict_res['TP']/(dict_res['TP'] + dict_res['FP'])
  try:
    spe = dict_res['FN']/(dict_res['TN'] + dict_res['FN'])
    return spe
  except ZeroDivisionError:
    return 'NA'



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
    
    # Guess the "mode":
    SAM_MODE = False
    if ext_infile == ".sam":
      SAM_MODE = True

    if SAM_MODE:
      print("SAM MODE !\n")

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

      # Set MAPQ cutoff:
      if num_cutoff == 'none':
        mapq_cutoff = 0
      else:
        mapq_cutoff = int(check.input_nb(num_cutoff, "-c cutoff (MAPQ here)"))

      
      
      print("Loading taxonomic Python module...")
      import src.shared_fun as shared
      tupl_sets_levels = shared.define(taxo_cutoff)      
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

      dict_gethered = {}
      dict_count = {"unmapped":0, "suppl_entries":0, "goods":0, "mapped":0}
      for alignment in input_samfile.fetch(until_eof=True):
        if alignment.is_unmapped:
            dict_count["unmapped"] += 1
            dict_stats['FN'] += 1 # Count unmapped = 'FN'
        else:
          if alignment.is_supplementary:
            dict_count["suppl_entries"] += 1
            dict_gethered[alignment.query_name].append(alignment)

          else:
            assert(alignment.query_name not in dict_gethered.keys())
            dict_gethered[alignment.query_name] = [alignment]
            
              
      input_samfile.close()
      print("SAM PARSING TIME:", str(t.time() - START_SAM_PARSING))


      # Handling supplementaries:
      print("\nHandling supplementary alignments...")
      # # Check if all lists have length >= 2:
      # assert(all(map(lambda a_list:len(a_list) > 1, dict_gethered.values())))

      list_res_suppl_handling = []
      list_MAPQ = []
      for query_name in dict_gethered:
        dict_count["mapped"] += 1
        align_list = dict_gethered[query_name]

        if len(align_list) > 1: # supplementary allignment
          res_suppl_handling = handle_suppl(align_list, dict_seqid2taxid)
          list_res_suppl_handling.append(res_suppl_handling)

        else: # Normal linear case
          align_obj = align_list[0]
          list_MAPQ.append(align_obj.mapping_quality)
          if align_obj.is_secondary: # We skip secondary
            pass
          elif align_obj.mapping_quality >= mapq_cutoff: #and align_obj.template_length > cutoff:
            # ref_name = alignment.reference_name
            # ref_name = alignment.reference_name.split()[0]
            #res_eval = eval_align(align_obj.reference_name, dict_seqid2taxid, 
            #                      taxo_cutoff)
            taxid = dict_seqid2taxid[align_obj.reference_name]
            res_eval = shared.in_zymo(taxid, taxo_cutoff, tupl_sets_levels)
            if res_eval == "TN":
              print("\n\nSALUT")
              print(shared.in_zymo(taxid, taxo_cutoff, tupl_sets_levels))
              print(align_obj.reference_name);sys.exit()
            # if res_eval == "TN": print("TNNNNEEG:", align_obj.reference_name, dict_seqid2taxid[align_obj.reference_name])
            dict_stats[res_eval] += 1
            dict_count["goods"] += 1
          else:
            dict_stats['FN'] += 1


      right_xlim = max(list_MAPQ)
      # assert(all(map(lambda x: x<=left_xlim, list_MAPQ)))
      plt.hist(list_MAPQ, bins=int(256/1), log=True)
      plt.xlim((-1, right_xlim))
      plt.show()
      
      # Print results of suppl handling:
      print(dict_count)
      tot_reads = dict_count["mapped"] + dict_count["unmapped"]
      print("% UNMAPPED:", dict_count["unmapped"]/tot_reads*100)
      print()

      nb_evaluated_suppl = len(list_res_suppl_handling)
      print("TOTAL EVALUATED (grouped) SUPPL:", nb_evaluated_suppl)
      mergeable_suppl = [elem for elem in list_res_suppl_handling if elem]
      print("MERGEABLE SUPPL:", len(mergeable_suppl))    
      print("REST (not mergeable):", nb_evaluated_suppl - len(mergeable_suppl))
      print()

      total_suppl_entries = dict_count["suppl_entries"] + nb_evaluated_suppl
      sum_stats, sum_counts = sum(dict_stats.values()),sum(dict_count.values())
      assert(sum_stats == tot_reads - nb_evaluated_suppl)
      print(sorted(dict_stats.items()))
      #print("TOT ALIGNMENTS:", sum_counts, " | ", "SUM STATS:", 
      #      sum_stats)
      
      dict_to_convert = dict_stats


    else:
      df_hits = pd.read_csv(to_infile, sep='\t', index_col=0, header=0)
      dict_stats_propag = {'TN':0, 'FN':0, 'TP':0, 'FP':0}

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
              
              dict_stats_propag[res] += df_hits.loc[clust, "nb_memb"]
              # if res == "FP":
              if False:
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
                  # print(dict_str[res])
                  dict_stats[res] += 1    
          
          
          else: # Not assigned because of the e_val cutoff
              print("Not assigned because of the e_val cutoff")
              dict_stats['FN'] += 1     
      
      print()
      print("AVANT PROPAG:", dict_stats)
      print("APRES PROPAG:", dict_stats_propag)
      assert(sum(dict_stats.values()) == len(rows_clust))
      assert(sum(dict_stats_propag.values()) == sum(df_hits["nb_memb"]))
      dict_to_convert = dict_stats_propag



    # COMPUTE METRICS:
    y_true, y_pred = dict_stats_to_vectors(dict_to_convert)
    sensitivity = round(skm.recall_score(y_true, y_pred), 4)
    NPV = 1 - calc_FOR(dict_to_convert)
    precision = round(skm.precision_score(y_true, y_pred), 4)

    print("SENSITIVITY:", sensitivity, " | ", "FNR:", 1 - sensitivity)
    print("FOR:", 1 - NPV, " | ", "NPV:", NPV)
    print("SPECIFICITY:", calc_specificity(dict_to_convert))
    print("PRECISION:", precision, " | ", "FDR:", 1 - precision)
    print("F1-SCORE:", round(skm.f1_score(y_true, y_pred), 4))
    # print(skm.classification_report(y_true, y_pred))
    print("MATTHEWS:", round(skm.matthews_corrcoef(y_true, y_pred), 4))
    print("ACCURACY:", round(skm.accuracy_score(y_true, y_pred), 4))
    print()
