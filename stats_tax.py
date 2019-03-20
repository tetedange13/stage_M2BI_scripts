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
            "AS":align_obj.get_tag("AS"),
            "is_second":align_obj.is_secondary}
            # "read_len":align_obj.infer_query_length()/align_obj.infer_read_length(),
            # "align_len":align_obj.infer_query_length()}


# def in_zymo_old(taxo_name, tupl_sets):
#     """
#     Given the rank of a BLAST hit and a sublist of taxonomic levels, deducted
#     from the cutoff (for the taxonomic level),  returns a string in 
#     ('TP', 'FP', 'TN', 'FN').
#     The name of the species is also requiered, to deal with 'strain' cases
#     """
#     # print(taxfoo.get_taxid_name(current_taxid), end='\t')
#     # current_lineage = taxfoo.get_lineage_as_dict(current_taxid)
#     # if taxonomic_cutoff not in current_lineage.keys():
#     #     return 'FP' # Taxo cutoff not in keys, so False Positive ?

#     # current_taxo_level = current_lineage[taxonomic_cutoff]
#     # assert(current_taxo_level)

#     set_levels_prok, set_levels_euk = tupl_sets
#     if taxo_name in set_levels_prok: # TP
#         return 'TP'
#     # elif taxo_name in set_levels_euk: # TN
#         # return 'TN'
#         # return 'FP' # Eukaryota from the Zymo treated as FP too 
#     else: # False Positive
#         return 'FP'


def str_from_res_eval(tupl_res_eval):
    """
    Generate a string (ready to write) from an evaluation result
    """
    readID, type_res, lineage, mapq = tupl_res_eval
    to_return = [readID, type_res]
    lineage_to_write = lineage 

    if type_res != 'ratio':
        lineage_to_write = ';'.join(str(taxid) for taxid in lineage)

    return ",".join(to_return + [lineage_to_write, str(mapq)])


def in_zymo(str_list_taxids, tupl_sets, taxonomic_cutoff):
    """
    """
    set_levels_prok, set_levels_euk = tupl_sets

    found = False
    for taxid in str_list_taxids.split(';'):
        if shared.taxfoo.get_taxid_rank(int(taxid)) == taxonomic_cutoff:
            found = True
            break

    taxo_name = shared.taxfoo.get_taxid_name(int(taxid))

    if not found:
        return 'FPnotInKey'
    else:
        if taxo_name in set_levels_prok:
            return (taxo_name, "TP")
    return (taxo_name, "FP")


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
    # Specificity or true negative rate: TNR = TN/(TN+FP) 
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


def compute_metrics(dict_stats_to_convert, at_taxa_level):
    """
    """
    import sklearn.metrics as skm
    y_true, y_pred = dict_stats_to_vectors(dict_stats_to_convert)

    precision = round(skm.precision_score(y_true, y_pred), 4) # PPV = TP/(TP+FP) (positive predictive value or precision)
    print("PRECISION:", precision, " | ", "FDR:", 1 - precision) 

    if not at_taxa_level:
        sensitivity = round(skm.recall_score(y_true, y_pred), 4)
    
        print("SENSITIVITY:", sensitivity, " | ", "FNR:", 1 - sensitivity) # TPR = TP/(TP+FN) (or sensitivity, hit rate, recall) 
        print("F1-SCORE:", round(skm.f1_score(y_true, y_pred), 4))
        print("MATTHEWS:", round(skm.matthews_corrcoef(y_true, y_pred), 4))
        print("ACCURACY:", round(skm.accuracy_score(y_true, y_pred), 4))
        print("SPECIFICITY:", calc_specificity(dict_stats_to_convert))
        # NPV = 1 - calc_FOR(dict_stats_to_convert)
        # print("FOR:", 1 - NPV, " | ", "NPV:", NPV)
    print()




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

        to_out_file = "salut.csv"
        if not osp.isfile(to_out_file):
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
            dict_count = {"mapped":0}
            list_suppl = []
            list_unmapped = []
            # dict_mapq = {} # Needed to filter out alignments with low MAPQ

            for idx, alignment in enumerate(input_samfile.fetch(until_eof=True)):
                if (idx+1) % 500000 == 0:
                    print("500 000 SAM entries elapsed !")
                query_name = alignment.query_name
                assert(query_name) # Different from ""

                if alignment.is_unmapped:
                    list_unmapped.append(query_name)
                    dict_stats['FN'] += 1 # Count unmapped = 'FN'

                else:
                    if alignment.has_tag("SA"): # We skip supplementaries
                        list_suppl.append(query_name)
                        if not alignment.is_supplementary:
                            dict_count["mapped"] += 1
                    else:
                        dict_align = alignment_to_dict(alignment)
                        if query_name in dict_gethered.keys():
                            dict_gethered[query_name].append(dict_align)
                        else:
                            # dict_mapq[query_name] = alignment.mapping_quality
                            # assert(query_name not in dict_gethered.keys())
                            dict_gethered[query_name] = [dict_align]
                  
                    
            input_samfile.close()
            print("SAM PARSING TIME:", str(t.time() - START_SAM_PARSING))

            assert(len(list_unmapped) == len(set(list_unmapped)))
            # dict_count['unmapped'] = len(list_unmapped)
            dict_count["mapped"] += len(dict_gethered.keys())
            dict_count["SA"] = len(list_suppl)
            dict_count["SA_uniq"] = len(set(list_suppl))
            del list_suppl

            # Removing supplementaries and low mapq:
            # fixed_mapq_cutoff = 2
            # all_query_name = list(dict_gethered.keys()) # List casting needed here
            # compt_suppl_SAM_entries, compt_nb_suppl = 0, 0

            # for query_name in all_query_name:
            #     len_align_list = len(dict_gethered[query_name])
            #     cond_1_a = len_align_list > 1
            #     # cond_2 = dict_gethered[query_name][0]["mapq"] < fixed_mapq_cutoff

            #     if (cond_1_a and dict_gethered[query_name][1]["is_suppl"]): #or cond_2:
            #         # if not cond_2:
            #         compt_suppl_SAM_entries += len_align_list
            #         compt_nb_suppl += 1
            #         print("TOTO")
            #         del dict_gethered[query_name]
            # print("NB_ENTRIES (SUPPL) REMOVED", compt_suppl_SAM_entries)
            # print("NB_SUPPL REMOVED", compt_nb_suppl)


            # Removing reads with lowest MAPQ:
            # CUTOFF_ON_MAPQ = 0 
            # # CUTOFF_ON_MAPQ = 0.05  # 5%

            # print("\nRemoving the", str(int(CUTOFF_ON_MAPQ*100)) + "% worst reads...")
            # sorted_by_mapq = sorted(dict_mapq.items(), key=lambda x: x[1])
            # inital_len_dict_gethered = len(sorted_by_mapq)
            # nb_to_remove = int(inital_len_dict_gethered*CUTOFF_ON_MAPQ)
            # print("(represents theorically:", nb_to_remove, "reads)")

            # i, nb_removed = 0, 0
            # while nb_removed < nb_to_remove and i < inital_len_dict_gethered:
            #     current_query_name = sorted_by_mapq[i][0]
            #     if len(dict_gethered[current_query_name]) == 1: # NOT suppl
            #         del dict_gethered[current_query_name]
            #         nb_removed += 1
            #     i += 1

            # if nb_removed < nb_to_remove:
            #     print("WARNING: Couldn't remove enough alignments...")
            #     print("(removed only", nb_removed)
            # assert(inital_len_dict_gethered == len(dict_gethered) + nb_removed)

            # if nb_removed != nb_to_remove:
            #     print("Number of alignments actually removed", nb_removed)
            # else:
            #     print("Well removed the predicted number of alignments !")
            print()


            # Handling supplementaries:
            TOTO = t.time()
            # CUTOFF_ON_RATIO = 0.9
            CUTOFF_ON_RATIO = 0
            

            print("Handling supplementary alignments...")
            print("Cutoff on alignment length of:", CUTOFF_ON_RATIO)
            # # Check if all lists have length >= 2:
            # assert(all(map(lambda a_list:len(a_list) > 1, dict_gethered.values())))
            # ref_name = alignment.reference_name
            # ref_name = alignment.reference_name.split()[0]
            partial_func = partial(pll.SAM_taxo_classif, 
                                   conv_seqid2taxid=dict_seqid2taxid,
                                   taxonomic_cutoff=taxo_cutoff, 
                                   tupl_sets=tupl_sets_levels,
                                   cutoff_ratio=CUTOFF_ON_RATIO)

            # Parallel version:
            my_pool = mp.Pool(15)
            results = my_pool.map(partial_func, dict_gethered.items())
            my_pool.close()
            # Serial version:
            # results = map(partial_func, dict_gethered.items())
    
            # Write outfile:
            with open(to_out_file, 'w') as salu_file:
                # Write header:
                salu_file.write(",type_align,lineage,mapq\n")

                # Write mapped reads:
                for res_eval in results:
                    salu_file.write(str_from_res_eval(res_eval) + '\n')

                # Write unmapped reads:
                for unmapped_read in list_unmapped:
                    salu_file.write(unmapped_read + ',unmapped,no\n')


            print("SUPPL HANDLING TIME:", t.time() - TOTO)
        # for res_eval in results:
        #     print(res_eval)
        #     readID, type_res = res_eval[0:2]

        #     if type_res == 'second' or type_res == 'linear':
        #         list_identified_sp.append(res_eval[2])
        #         if type_res == 'second':
        #             dict_count["secondaries"] += 1
        #     elif type_res == 'FP':
        #         dict_stats[type_res] += 1
        #         compt_FP += 1
        #         list_MAPQ.append(res_eval[1])
        #     else: # type_res 'ratio':
        #         ratio = res_eval[1]
                    
        else:
            print("FOUND CSV FILE !")

        # list_MAPQ = [] # Contain only MAPQ of align that pass all filters
        # list_identified_sp = []
        # compt_FP = 0

        print("Loading CSV file...")
        my_csv = pd.read_csv(to_out_file, header=0, index_col=0)
        print("CSV loaded !")
        dict_count = {"unmapped":0}

        for readID in my_csv.index:
            type_align = my_csv.loc[readID, "type_align"]
            if type_align == "unmapped":
                dict_count['unmapped'] += 1
            elif type_align == 'ratio':
                pass
            else: # normal or secondary
                species, res = in_zymo(my_csv.loc[readID, "lineage"], 
                                              tupl_sets_levels, taxo_cutoff)
                print(species, res)


        sys.exit()
        # print("FP FROM NOT IN KEYS:", compt_FP)

        set_seen_prok = set()
        dict_species2res = {} # To access evaluation results of a given species 

        for species in list_identified_sp:
            res_eval = in_zymo(species, tupl_sets_levels)
            dict_stats[res_eval] += 1
            dict_species2res[species] = res_eval

            if res_eval == 'TP': # If it is a proK from th Zymo
                set_seen_prok.add(species)
        

         # AT THE TAXA LEVEL:
        # set_prok_Zymo = tupl_sets_levels[0]
        recall_at_taxa_level = len(set_seen_prok)/len(tupl_sets_levels[0])
        # assert(all(map(lambda dict_val: dict_val != 'FN', 
        #                dict_species2res.values())))

        print("RESULTS AT THE TAXA LEVELS:")
        TP_list = [sp for sp in dict_species2res if dict_species2res[sp] == 'TP']
        print("NB_TP", len(TP_list), " | NB_FP", len(dict_species2res))
        dict_stats_sp_level = {'TP':len(TP_list),
                               'FP':len(dict_species2res.keys())}
        compute_metrics(dict_stats_sp_level, True)
        print("SENSITIVITY:", recall_at_taxa_level, " | ", "FNR:",
              1 - recall_at_taxa_level)


        if list_MAPQ:
            # Display distribution of MAPQ:
            right_xlim = max(list_MAPQ)
            # assert(all(map(lambda x: x<=left_xlim, list_MAPQ)))
            plt.hist(list_MAPQ, bins=int(256/1), log=True)
            plt.xlim((-1, right_xlim))
            plt.show()
        
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
    compute_metrics(dict_to_convert, False)
