#!/usr/bin/env python3


"""
Mapping statistics computation

Usage:
  main.py (-i <inFile>) (-l <taxoCut>)
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inFile=input_file     input file
  -l --taxoCut=taxo_cutoff   cutoff for the taxonomic level [default: none]
"""


import sys, os
import os.path as osp
import pandas as pd
import time as t
import subprocess as sub
import pysam as pys
import multiprocessing as mp
import sklearn.metrics as skm
import matplotlib.pyplot as plt
from functools import partial
from docopt import docopt
import src.remote as remote
import src.check_args as check
import src.parallelized as pll


def soft_clipping_from_cigar(list_cigar_tupl):
    """
    """
    start_soft_clip, end_soft_clip = list_cigar_tupl[0], list_cigar_tupl[-1]

    nb_start_soft_clip, nb_end_soft_clip = 0, 0
    if start_soft_clip[0] == 4:
        nb_start_soft_clip = start_soft_clip[1]
    if end_soft_clip[0] == 4:
        nb_end_soft_clip = end_soft_clip[1]

    return (nb_start_soft_clip, nb_end_soft_clip)



def alignment_to_dict(align_obj):
    """
    Takes an alignment instance and return a dict containing the needed
    attributes (only)

    With infer_read_length() method, hard-clipped bases are included in the 
    counting
    infer_query_length() method does NOT include them
    """
    ratio_len = align_obj.query_length/align_obj.infer_read_length()
    assert(ratio_len <= 1)

    return {"mapq":align_obj.mapping_quality,
            "ref_name": align_obj.reference_name,
            "ratio_len":ratio_len,
            "is_suppl":align_obj.is_supplementary,
            "is_second":align_obj.is_secondary}
            # "read_len":align_obj.infer_query_length()/align_obj.infer_read_length(),
            # "align_len":align_obj.infer_query_length()}


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
    # Sensitivity, hit rate, recall, or true positive rate: TPR = TP/(TP+FN)
    # Specificity or true negative rate: TNR = TN/(TN+FP) 
    # Precision or positive predictive value: PPV = TP/(TP+FP)
    # Negative predictive value: NPV = TN/(TN+FN)
    # Fall out or false positive rate: FPR = FP/(FP+TN)
    # False negative rate: FNR = FN/(TP+FN)
    # False discovery rate: FDR = FP/(TP+FP)
    # Overall accuracy: ACC = (TP+TN)/(TP+FP+FN+TN)
    try:
        spe = dict_res['TN']/(dict_res['TN'] + dict_res['FP'])
        return spe
    except ZeroDivisionError:
        return 'NA'


def calc_FOR(dict_res):
    # Negative predictive value: NPV = TN/(TN+FN)
    # FOR = 1 - NPV

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
    to_infile, infile_base, _, _, ext_infile = check.infile(ARGS["--inFile"], 
                                                         ["csv", "tsv", "sam"])
    taxo_cutoff = ARGS["--taxoCut"]

    # Common variables:
    to_apps = "/home/sheldon/Applications/"
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    # taxo_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 
    #                'genus', 'species', 'subspecies'] + ["strain"] 
    dict_stats = {'TN':0, 'FN':0, 'TP':0, 'FP':0}

    print("Loading taxonomic Python module...")
    import src.shared_fun as shared
    tupl_sets_levels = shared.generate_sets_zymo(taxo_cutoff)      
    print("Taxonomic Python module loaded !\n")
    
    # Guess the "mode":
    SAM_MODE = False
    if ext_infile == ".sam":
        SAM_MODE = True

    if SAM_MODE:
        print("SAM MODE !")

        # Guess the database used:
        root_to_seqid2taxid = to_dbs + "Centri_idxes/"
        if "toZymo" in infile_base:
            print("DB GUESSED: Zymo")
            root_to_seqid2taxid += "Zymo"
        elif "toRrn" in infile_base:
            print("DB GUESSED: rrn")
            root_to_seqid2taxid += "rrn"
        elif "toSilva" in infile_base:
            print("DB GUESSED: SILVA")
            root_to_seqid2taxid += "SILVA"
        else:
            print("Unkown database !\n")
            sys.exit(2)
        to_seqid2taxid = root_to_seqid2taxid + "/seqid2taxid"
        print()

        # To make correspond operon number and taxid:
        dict_seqid2taxid = {}
        with open(to_seqid2taxid, 'r') as seqid2taxid_file:
            for line in seqid2taxid_file:
                splitted_line = line.rstrip('\n').split('\t')
                dict_seqid2taxid[splitted_line[0]] = int(splitted_line[1])

        # Start SAM file parsing:
        print("Extracting information from SAM file...")
        START_SAM_PARSING = t.time()
        input_samfile = pys.AlignmentFile(to_infile, "r")

        dict_gethered = {}
        dict_count = {"unmapped":0, "suppl_entries":0, "goods":0, "mapped":0}
        dict_mapq = {} # Needed to filter out alignments with low MAPQ

        for idx, alignment in enumerate(input_samfile.fetch(until_eof=True)):
            if (idx+1) % 500000 == 0: print("500 000 SAM entries elapsed !")

            if alignment.is_unmapped:
                dict_count["unmapped"] += 1
                dict_stats['FN'] += 1 # Count unmapped = 'FN'

            else:
                # print(soft_clipping_from_cigar(alignment.cigartuples))
                dict_align = alignment_to_dict(alignment)
                query_name = alignment.query_name
                if alignment.is_supplementary or alignment.is_secondary:
                    # dict_count["suppl_entries"] += 1
                    dict_gethered[query_name].append(dict_align)
                else:
                    dict_mapq[query_name] = alignment.mapping_quality
                    assert(query_name not in dict_gethered.keys())
                    dict_gethered[query_name] = [dict_align]
              
                
        input_samfile.close()
        print("SAM PARSING TIME:", str(t.time() - START_SAM_PARSING))
        dict_count["mapped"] = len(dict_gethered.keys())

        # Removing reads with lowest MAPQ:
        # CUTOFF_ON_MAPQ = 0.1  # 10%
        CUTOFF_ON_MAPQ = 0.05  # 5%

        print("\nRemoving the", str(int(CUTOFF_ON_MAPQ*100)) + "% worst reads...")
        sorted_by_mapq = sorted(dict_mapq.items(), key=lambda x: x[1])
        inital_len_dict_gethered = len(sorted_by_mapq)
        nb_to_remove = int(inital_len_dict_gethered*CUTOFF_ON_MAPQ)
        print("(represents theorically:", nb_to_remove, "reads)")

        i, nb_removed = 0, 0
        while nb_removed < nb_to_remove and i < inital_len_dict_gethered:
            current_query_name = sorted_by_mapq[i][0]
            if len(dict_gethered[current_query_name]) == 1: # NOT suppl
                del dict_gethered[current_query_name]
                nb_removed += 1
            i += 1

        if nb_removed < nb_to_remove:
            print("WARNING: Couldn't remove enough alignments...")
            print("(removed only", nb_removed)
        assert(inital_len_dict_gethered == len(dict_gethered) + nb_removed)

        if nb_removed != nb_to_remove:
            print("Number of alignments actually removed", nb_removed)
        else:
            print("Well removed the predicted number of alignments !")
        print()
        # dict_stats['FN'] += nb_removed


        # Handling supplementaries:
        TOTO = t.time()
        CUTOFF_ON_RATIO = 0.9
        # CUTOFF_ON_RATIO = 0
        

        print("\nHandling supplementary alignments...")
        print("Cutoff on alignment length of:", CUTOFF_ON_RATIO)
        # # Check if all lists have length >= 2:
        # assert(all(map(lambda a_list:len(a_list) > 1, dict_gethered.values())))
        # ref_name = alignment.reference_name
        # ref_name = alignment.reference_name.split()[0]

        # Parallel version:
        partial_func = partial(pll.SAM_taxo_classif, 
                               conv_seqid2taxid=dict_seqid2taxid,
                               taxonomic_cutoff=taxo_cutoff, 
                               tupl_sets=tupl_sets_levels,
                               cutoff_ratio=CUTOFF_ON_RATIO)
        my_pool = mp.Pool(10)
        results = my_pool.map(partial_func, dict_gethered.values())
        my_pool.close()

        # Serial version:
        # results = []
        # for query_name in dict_gethered:
        #   results.append(pll.SAM_taxo_classif(dict_gethered[query_name], 
        #                  dict_seqid2taxid, taxo_cutoff, tupl_sets_levels, 
        #                  CUTOFF_ON_RATIO))

        # list_res_suppl_handling = []
        list_MAPQ = [] # Contain only MAPQ of align that pass all filters
        for res_eval in results:
            first_elem = res_eval[0]
            if first_elem == 'suppl':
                pass
                # list_res_suppl_handling.append(res_eval[1])
            else:
                # if res_eval[0] == 'FN':
                #     print("VAL RATIO:", res_eval[1]);sys.exit()
                if first_elem != 'ratio':
                    dict_stats[res_eval[0]] += 1
                    list_MAPQ.append(res_eval[1])
        print("SUPPL HANDLING TIME:", t.time()-TOTO)

        # Display distribution of MAPQ:
        # right_xlim = max(list_MAPQ)
        # # assert(all(map(lambda x: x<=left_xlim, list_MAPQ)))
        # plt.hist(list_MAPQ, bins=int(256/1), log=True)
        # plt.xlim((-1, right_xlim))
        # plt.show()
        
        # Print general counting results:
        print(dict_count)
        tot_reads = dict_count["mapped"] + dict_count["unmapped"]
        print("TOT READS:", tot_reads)
        print("% UNMAPPED:", dict_count["unmapped"]/tot_reads*100)
        print()

        # Print results of suppl handling:
        # nb_evaluated_suppl = len(list_res_suppl_handling)
        # total_suppl_entries = dict_count["suppl_entries"] + nb_evaluated_suppl
        # print("TOTAL EVALUATED (grouped) SUPPL:", nb_evaluated_suppl)
        # mergeable_suppl = [elem for elem in list_res_suppl_handling if elem]
        # print("MERGEABLE SUPPL:", len(mergeable_suppl))    
        # print("REST (not mergeable):", nb_evaluated_suppl - len(mergeable_suppl))
        # print()

        
        sum_stats, sum_counts = sum(dict_stats.values()),sum(dict_count.values())
        # assert(sum_stats == tot_reads - nb_evaluated_suppl)
        print(sorted(dict_stats.items()))
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
                res = shared.in_zymo(df_hits.loc[clust, "taxid"], taxo_cutoff, 
                                   tupl_sets_levels)
                # res = in_zymo(df_hits.loc[clust, "topHit"], 
                #               df_hits.loc[clust, "rank"], 
                #               taxo_cutoff)
              
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
                                  " | ", "RANK:", 
                                  df_hits.loc[alt_index, "rank"]) 
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
    print("SPECIFICITY:", calc_specificity(dict_to_convert))
    print("PRECISION:", precision, " | ", "FDR:", 1 - precision)
    print("F1-SCORE:", round(skm.f1_score(y_true, y_pred), 4))
    # print(skm.classification_report(y_true, y_pred))
    print("FOR:", 1 - NPV, " | ", "NPV:", NPV)
    print("MATTHEWS:", round(skm.matthews_corrcoef(y_true, y_pred), 4))
    print("ACCURACY:", round(skm.accuracy_score(y_true, y_pred), 4))
    print()
