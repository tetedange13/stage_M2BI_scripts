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


def alignment_to_dict(align_obj):
    """
    Takes an alignment instance and return a dict containing the needed
    attributes (only)

    With infer_read_length() method, hard-clipped bases are included in the 
    counting
    infer_query_length() method does NOT include them
    """
    ratio_len = align_obj.infer_query_length()/align_obj.infer_read_length()
    assert(ratio_len <= 1)

    to_return = {"mapq":align_obj.mapping_quality,
                 "ref_name": align_obj.reference_name,
                 "len_align":align_obj.infer_query_length(),
                 "ratio_len":ratio_len,
                 "de":round(align_obj.get_tag("de"), 4),
                 "is_suppl":align_obj.is_supplementary,
                 "has_SA": align_obj.has_tag("SA"),
                 "AS":align_obj.get_tag("AS"),
                 "is_second":align_obj.is_secondary}

    return to_return


def plot_thin_hist(list_values, title_arg="", y_log=True, xlims=(0.15, 0.3)):
    """
    Draw a thin histogram from a list of values
    """
    fig = plt.figure()
    axis = plt.subplot(111)
    plt.title(title_arg)

    max_val = max(list_values)
    min_val = min(list_values)


    plt.hist(list_values, bins=int(256/1), log=y_log)
    plt.xlim(xlims)
    plt.show()


def str_from_res_eval(tupl_res_eval):
    """
    Generate a string (ready to write) from an evaluation result
    """
    readID, type_res, lineage = tupl_res_eval[0:3]
    rest = tupl_res_eval[3: ]
    to_return = [readID, type_res]
     
    if type_res == 'ratio' or type_res == 'pb_lca':
        lineage_to_write = lineage
    else:
        lineage_to_write = lineage
        # lineage_to_write = ';'.join(str(taxid) for taxid in lineage)

    return ",".join(to_return + [lineage_to_write] + [str(x) for x in rest])


def calc_taxo_shift(arg_taxid, taxonomic_cutoff):
    """
    Calculate the 'taxonomic shift', i.e. the number of taxonomic levels of
    difference, between a taxonomic cutoff and the rank of a taxid
    """
    idx_cutoff = shared.want_taxo.index(taxonomic_cutoff) + 1
    sublist_taxo = shared.want_taxo[0:idx_cutoff]
    shift = len(shared.taxfoo.get_lineage(arg_taxid)) - len(sublist_taxo)

    print("BONJ", shared.taxfoo.get_taxid_name(arg_taxid), shift)
    return shift
    # if shift <= 0:
    #     pass
    
    # print(len(shared.taxfoo.))


def old_lca_last_try(str_list_taxids):
    """
    Try to produce a LCA from a set of taxid, by eliminating
    """
    # old_set_taxids = map(int, str_list_taxids.split('&'))
    old_set_taxids = list(map(int, str_list_taxids.split('&')))
    new_set_taxids = set()
    list_ranks = []
    problems = 'aucun'
    for taxid in old_set_taxids:
        taxid_name = shared.taxfoo.get_taxid_name(taxid)
        if ('metagenome' in taxid_name or 'uncultured' in taxid_name or 
            'unidentified' in taxid_name or 'phage' in taxid_name.lower() or
            taxid_name == 'synthetic construct' or 
            'virus' in taxid_name.lower()):
            pass
            # print("PB:", taxid_name, taxid)
        else:
            new_set_taxids.add(taxid)

    if not new_set_taxids:
        return 'only_trashes'
    return shared.taxfoo.find_lca(new_set_taxids)


def lca_last_try(old_set_taxids):
    """
    Try to produce a LCA from a set of taxid, by eliminating
    """
    # old_set_taxids = list(map(int, str_list_taxids.split('&')))
    new_set_taxids = set()
    list_ranks = []
    for taxid in old_set_taxids:
        taxid_name = shared.taxfoo.get_taxid_name(taxid)
        if ('metagenome' in taxid_name or 'uncultured' in taxid_name or 
            'unidentified' in taxid_name or 'phage' in taxid_name.lower() or
            taxid_name == 'synthetic construct' or 
            'virus' in taxid_name.lower()):
            pass
            # print("PB:", taxid_name, taxid)
        else:
            new_set_taxids.add(taxid)

    if not new_set_taxids:
        return 'only_trashes'
    return shared.taxfoo.find_lca(new_set_taxids)


def handle_second(str_list_taxids):
    """
    Take an alignment that has been found to be secondary and try to determine
    its taxonomy (i.e. 1 unique taxid)    
    """
    set_taxid_target = set(map(int, str_list_taxids.split('s')))
    lca = shared.taxfoo.find_lca(set_taxid_target)

    if lca == 1: # If LCA searching fails, LCA==1
        lca_attempt = lca_last_try(set_taxid_target)
        if lca_attempt == 'only_trashes':
            return (lca_attempt, )
        elif lca_attempt == 1:
            return ('unsolved_lca_pb', )
        else:
            return ('lca_last_try', lca_attempt)

    return ('lca', lca)


def in_zymo(taxo_taxid, tupl_sets, taxonomic_cutoff):
    """
    Given the taxid of a read (lowest one in the taxo), determine if the 
    organism belongs to the Zymo mock comm, at a given taxonomic cutoff 
    """
    set_levels_prok, set_levels_euk = tupl_sets
    lineage = shared.taxfoo.get_lineage_as_dict(taxo_taxid)
    found = False
    taxo_levels = lineage.keys()

    if taxonomic_cutoff not in taxo_levels:
        return ('notDeterminable', 'FP')
    else:
        taxo_name = lineage[taxonomic_cutoff]
        if taxo_name in set_levels_prok:
            return (taxo_name, "TP")
        return (taxo_name, "FP")


def final_eval(csv_index_val, two_col_from_csv, sets_levels, taxonomic_cutoff):
    """
    """
    type_align = two_col_from_csv.loc[csv_index_val, "type_align"]

    # if type_align in ('ratio', 'unmapped'):
    #     pass
    # else: # normal or secondary
    lineage_val = two_col_from_csv.loc[csv_index_val, "lineage"]
    
    if type_align != 'pb_lca':
        taxid_to_eval = int(lineage_val.strip(';'))
    else:
        lca_attempt = lca_last_try(lineage_val)
        if lca_attempt == 'only_trashes':
            return (csv_index_val, 'only_trashes', 'only_trashes', 'FP')
        elif lca_attempt == 1:
            return (csv_index_val, 'unsolved_lca_pb', 'unsolved_lca_pb', 'FP')
        taxid_to_eval = lca_attempt
        # calc_taxo_shift(lca_attempt, taxonomic_cutoff)
        # else:
        #     lineage_to_eval = shared.taxfoo.get_lineage_as_taxids(lca_attempt)

    taxo_name, res_evaluation = in_zymo(taxid_to_eval, sets_levels, 
                                        taxonomic_cutoff)
    return (csv_index_val, taxid_to_eval, taxo_name, res_evaluation)


def eval_taxo(csv_index_val, two_col_from_csv, sets_levels, taxonomic_cutoff):
    """
    """
    lineage_val = two_col_from_csv.loc[csv_index_val, "lineage"]
    type_align = two_col_from_csv.loc[csv_index_val, "type_align"]
    if lineage_val.startswith(';'): # normal or second_uniq
        taxid_to_eval = lineage_val.strip(';')
    else: # Secondary with more than 1 unique taxid
        res_second_handling = handle_second(lineage_val)
        remark_eval = res_second_handling[0]
        if len(res_second_handling) == 1: # Problem (only trashes or unsolved)
            return (csv_index_val, remark_eval, 'FP', remark_eval)
        else:
            taxid_to_eval = res_second_handling[1]

    taxo_name, classif = in_zymo(taxid_to_eval, sets_levels, taxonomic_cutoff)
    # return (csv_index_val, taxid_to_eval, taxo_name, classif)

    if lineage_val.startswith(';'):
        return (csv_index_val, taxo_name, classif, type_align)
    return (csv_index_val, taxo_name, classif, remark_eval)


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
        print("ACCURACY:", round(skm.accuracy_score(y_true, y_pred), 4)) # ACC = (TP+TN)/(TP+FP+FN+TN) (overall accuracy)
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

        to_out_file = infile_base + ".csv"
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
            list_suppl = []
            list_unmapped = []

            for idx, alignment in enumerate(input_samfile.fetch(until_eof=True)):
                if (idx+1) % 500000 == 0:
                    print("500 000 SAM entries elapsed !")
                query_name = alignment.query_name
                assert(query_name) # Different from ""

                if alignment.is_unmapped:
                    list_unmapped.append(query_name)
                    dict_stats['FN'] += 1 # Count unmapped = 'FN'

                else:
                    dict_align = alignment_to_dict(alignment)
                    if alignment.is_secondary or alignment.is_supplementary:
                        dict_gethered[query_name].append(dict_align)
                    else:
                        assert(query_name not in dict_gethered.keys())
                        dict_gethered[query_name] = [dict_align]
                  
                    
            input_samfile.close()
            print("SAM PARSING TIME:", str(t.time() - START_SAM_PARSING))
            print()
            assert(len(list_unmapped) == len(set(list_unmapped)))
            list_suppl_as_repr = []

            # print("SUPPL AS REPR:", len(list_suppl_as_repr))
            # plot_thin_hist(list_nb_second, 
            #                "Distrib of nb of secondaries with N=100")
            # sys.exit()
            # dict_count["SA"] = len(list_suppl)
            # dict_count["SA_uniq"] = len(set(list_suppl))
            # del list_suppl


            # Extract releval SAM info towards a CSV file:
            TOTO = t.time()
            print("Converting SAM into CSV...")
            partial_func = partial(pll.SAM_to_CSV, 
                                   conv_seqid2taxid=dict_seqid2taxid)
                                               # Parallel version:
            my_pool = mp.Pool(15)
            results = my_pool.map(partial_func, dict_gethered.items())
            my_pool.close()
            # Serial version (need list casting to have output):
            # results = list(map(partial_func, dict_gethered.items())) 

            # sys.exit()

            # Write outfile:
            with open(to_out_file, 'w') as out_file:
                # Write header:
                out_file.write(",type_align,lineage,nb_trashes,mapq,len_align," +
                               "ratio_len,de\n")

                # Write mapped reads:
                for res_eval in results:
                    out_file.write(str_from_res_eval(res_eval) + '\n')

                # Write unmapped reads:
                for unmapped_read in list_unmapped:
                    out_file.write(unmapped_read + ',unmapped,no\n')

            print("CSV CONVERSION TIME:", t.time() - TOTO)
            sys.exit()

                    
        else:
            print("FOUND CSV FILE !")


        TIME_CSV_TREATMENT = t.time()
        print("Loading CSV file...")
        my_csv = pd.read_csv(to_out_file, header=0, index_col=0)
        print("CSV loaded !")

        # list_lineage_second = my_csv[my_csv['type_align'] == 'second_plural']['lineage']
        # list_second_len = list(map(lambda val: len(val.split('s')), list_lineage_second))
        # list_second_len = list(map(lambda val: len(set(val.split('s'))), list_lineage_second))
        # plot_thin_hist(list_second_len, "Distrib len_second", False, (-1, 12))
        # sys.exit()

        with_lineage = ((my_csv['type_align'] != 'unmapped') & 
                (my_csv['type_align'] != 'only_suppl'))
        my_csv_to_pll = my_csv[['lineage', 'type_align']][with_lineage]
        partial_eval = partial(eval_taxo, two_col_from_csv=my_csv_to_pll,
                                          sets_levels=tupl_sets_levels,
                                          taxonomic_cutoff=taxo_cutoff)

        # Parallel version:
        eval_pool = mp.Pool(15)
        results_eval = eval_pool.map(partial_eval, my_csv_to_pll.index)
        eval_pool.close()
        del my_csv_to_pll
        print("PLL PROCESS FINISHED !")
        
        # Compute counts:
        nb_unmapped = sum(my_csv['type_align'] == 'unmapped')
        nb_second_plur = sum(my_csv['type_align'] == 'second_plural')
        nb_second_uniq = sum(my_csv['type_align'] == 'second_uniq')
        dict_count = {"unmapped":nb_unmapped,
                      "only_suppl":sum(my_csv['type_align'] == 'only_suppl'),
                      "normal":sum(my_csv['type_align'] == 'normal'),
                      'second_plural':nb_second_plur,
                      'second_uniq':nb_second_uniq,
                      'tot_second':nb_second_plur + nb_second_uniq,
                      "unsolved_lca_pb":0, 'FPnotInKey':0, "only_trashes":0}
        dict_count["tot_reads"] = (nb_unmapped + dict_count["tot_second"] +
                                   dict_count["normal"])
        dict_stats['FN'] += nb_unmapped

        set_seen_prok = set()
        dict_species2res = {} # To access evaluation results of a given species
        list_index_new_col, list_val_new_col = [], []
        problems = ('only_trashes', 'notDeterminable', 'unsolved_lca_pb')

        for tupl_res in results_eval:
            readID, species, res_eval, remark_evaluation = tupl_res

            if species not in problems:
                dict_species2res[species] = res_eval
            else:
                if species == 'unsolved_lca_pb':
                    dict_count["unsolved_lca_pb"] += 1
                elif species == 'notDeterminable':
                    dict_count['FPnotInKey'] += 1
                else:
                    dict_count['only_trashes'] += 1

                if res_eval == 'TP':
                    set_seen_prok.add(species)
                    
            list_index_new_col.append(readID)
            list_val_new_col.append((species, res_eval, remark_evaluation))

        tmp_df = pd.DataFrame({'species':[tupl[0] for tupl in list_val_new_col],
                               'res_eval':[tupl[1] for tupl in list_val_new_col],
                               'remark_eval':[tupl[2] for tupl in list_val_new_col]},
                              index=list_index_new_col)
        my_csv = my_csv.assign(species=tmp_df['species'], 
                               res_eval=tmp_df['res_eval'],
                               remark_eval=tmp_df['remark_eval'])
        del tmp_df
        print("TIME FOR CSV TREATMENT:", t.time() - TIME_CSV_TREATMENT)
        print()

        toto = my_csv[my_csv["species"] == 'notDeterminable']
        for val in toto.index:
            if toto.loc[val, "nb_trashes"] == 0:
                # print(toto.loc[val, "type_align"])
                # print(toto.loc[val, "lineage"])
                pass

        # list_MAPQ = my_csv['mapq'][my_csv['mapq'].notnull()]
        # plot_thin_hist(list_MAPQ, "Distrib MAPQ", True, (-1, 60))        

        # FIXED_MAPQ_CUTOFF = 0
        # print(sum(map(lambda mapq_val: mapq_val>=FIXED_MAPQ_CUTOFF, 
        #               my_csv["mapq"])))
        # CUTOFF_ON_RATIO = 0.9
        # CUTOFF_ON_RATIO = 0
        # print("Cutoff on alignment length of:", CUTOFF_ON_RATIO)

        # cutoff_nb_trashes = 4
        cutoff_nb_trashes = max(my_csv['nb_trashes'])
        print("CUTOFF ON THE NB OF TRASHES:", cutoff_nb_trashes)
        print()

        is_FP = ((my_csv['res_eval'].notnull()) & 
                 (my_csv['res_eval'] == 'FP') &
                 (my_csv['nb_trashes'] <= cutoff_nb_trashes))
        is_TP = ((my_csv['res_eval'].notnull()) & 
                 (my_csv['res_eval'] == 'TP') &
                 (my_csv['nb_trashes'] <= cutoff_nb_trashes))
        print("FP STATS:")
        print(my_csv[is_FP]['nb_trashes'].value_counts().sort_index())
        # print(my_csv[is_FP & (my_csv['species'] == 'notDeterminable')]['nb_trashes'].value_counts())
        print()
        print("TP STATS:")
        # print(my_csv[is_TP]['nb_trashes'].value_counts().sort_index())
        # plt.figure()
        # axis = plt.subplot(111)
        # plt.title("SALUT")
        pd.DataFrame({'TP_nb_trashes':my_csv[is_TP]['nb_trashes'].value_counts(), 
                      'FP_nb_trashes':my_csv[is_FP]['nb_trashes'].value_counts()
                      }).plot.bar(title="Distrib nb_trashes between FP and TP",
                                  color=('red', 'green'))
        plt.show()
        # my_csv[is_FP]['nb_trashes'].value_counts().plot(kind='bar', 
        #                                                 colorbar=False,
        #                                                 )
        print()

        dict_stats['FP'] += sum(is_FP)
        dict_stats['FP'] += sum(my_csv['type_align'] == 'only_suppl')
        dict_stats['TP'] += sum(is_TP)
        print(dict_count)
        print(dict_stats)
        print("TOT READS EVALUATED:", sum(dict_stats.values()))
        print()

        # true_FP = ((my_csv['res_eval'] == 'FP') & 
        #            (my_csv['species'] != 'notDeterminable') &
        #            (my_csv['species'] != 'unsolved_lca_pb'))
        
        # print(my_csv['species'][false_pos].)
        # salut = shared.taxfoo.get_lineage_as_dict(list_taxids[-1])
        # print(my_csv[is_FP & (my_csv['species'] == 'notDeterminable')]['type_align'].value_counts())
        salut = my_csv[is_FP & 
                       (my_csv['species'] == 'notDeterminable') &
                       (my_csv['type_align'] == 'merged')]['lineage']
        for taxid_lineage in salut:
            # print(taxid_lineage.strip(';'))
            calc_taxo_shift(int(taxid_lineage.strip(';')), taxo_cutoff)
            # lca = lca_last_try(taxid_lineage)
            # if lca != 'only_trashes':
            #     print(lca)

                # if lca == 204429:
                #     print(shared.taxfoo.get_lineage_as_dict(lca))

            # calc_taxo_shift(val, taxo_cutoff)
            # calc_taxo_shift()
            #print(val)
        
        # plot_thin_hist(my_csv[is_TP]['mapq'], "Distrib MAPQ", True, (-1, 60))

        # print(my_csv[is_TP]['type_align'].value_counts())
        # print()
        # print(my_csv[is_FP & 
        #      (my_csv['species']=='notDeterminable') & 
        #      (my_csv['type_align'] == 'merged')]['mapq'].value_counts())
        # print(my_csv[is_FP]['species'].value_counts())
        # print(my_csv[false_pos].groupby('species').count())


        sys.exit()


        # list_len_align = []
        # plot_thin_hist(list_len_align, "Distrib len_align TP", True)
        


        # AT THE TAXA LEVEL:
        print("RESULTS AT THE TAXA LEVELS:")
        recall_at_taxa_level = len(set_seen_prok)/len(tupl_sets_levels[0])
        TP_list = [sp for sp in dict_species2res if dict_species2res[sp] == 'TP']
        print("NB_TP", len(TP_list), " | NB_FP", len(dict_species2res))
        dict_stats_sp_level = {'TP':len(TP_list),
                               'FP':len(dict_species2res.keys())}
        print("SENSITIVITY:", recall_at_taxa_level, " | ", "FNR:",
              1 - recall_at_taxa_level)
        compute_metrics(dict_stats_sp_level, True)
        


        # Print general counting results:
        tot_reads = dict_count['total_reads']
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
