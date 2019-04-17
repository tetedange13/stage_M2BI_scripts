#!/usr/bin/env python3


"""
Mapping statistics computation

Usage:
  stats_tax.py (-i <inFile>) (-l <taxoCut>)
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inFile=input_file     input file
  -l --taxoCut=taxo_cutoff   cutoff for the taxonomic level
"""


import sys, os
import os.path as osp
import time as t
import subprocess as sub
import pysam as pys
import multiprocessing as mp
import matplotlib.pyplot as plt
from functools import partial
from docopt import docopt
import src.check_args as check


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
                 "is_suppl":align_obj.is_supplementary,
                 "has_SA": align_obj.has_tag("SA"),
                 "AS":align_obj.get_tag("AS"),
                 "is_second":align_obj.is_secondary}
                 
    if align_obj.has_tag('de'):
        to_return['de'] = round(align_obj.get_tag("de"), 4)

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


def str_from_res_conv(dict_res_conv):
    """
    Generate a string (ready to write) from an evaluation result
    """
    to_extract = ['readID', 'type_align', 'lineage']
    # ,type_align,lineage,nb_trashes,mapq,len_align," +
                #                "ratio_len,de\n")
    rest = [(a_key, dict_res_conv[a_key]) for a_key in dict_res_conv 
                                          if a_key not in to_extract]
    to_return = [dict_res_conv['readID'], dict_res_conv['type_align']]
     
    if ( dict_res_conv['type_align'] == 'ratio' or 
         dict_res_conv['type_align'] == 'pb_lca' ):
        lineage_to_write = dict_res_conv['lineage']
    else:
        lineage_to_write = dict_res_conv['lineage']
        # lineage_to_write = ';'.join(str(taxid) for taxid in lineage)

    return (",".join(to_return + [lineage_to_write] + 
                     [str(x[1]) for x in rest]), 
            to_extract[1:] + list(map(lambda x: x[0], rest)))


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
    else:
        return precision



# MAIN:
if __name__ == "__main__":
    # ARGUMENTS:
    ARGS = docopt(__doc__, version='0.1')
    to_infile, infile_base, _, _, ext_infile = check.infile(ARGS["--inFile"], 
                                                ['csv', 'tsv', 'txt', 'sam'])

    # Common variables:
    to_apps = "/home/sheldon/Applications/"
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    dict_stats = {'TN':0, 'FN':0, 'TP':0, 'FP':0}

    print("Loading taxonomic Python module...")
    import src.parallelized as pll
    evaluate = pll.eval
    pd = evaluate.pd
    print("Taxonomic Python module loaded !\n")

    taxo_cutoff = check.acceptable_str(ARGS["--taxoCut"], evaluate.want_taxo)
    # Guess the database used:
    if "toZymo" in infile_base:
        guessed_db = 'zymo'
    elif "toRrn" in infile_base:
        guessed_db = 'rrn'
    elif "toSilva" in infile_base:
        guessed_db = 'silva'
    elif "toP_compressed" in infile_base:
        guessed_db = 'p_compressed'
    else:
        print("Unkown database !\n")
        sys.exit(2)
    print("DB GUESSED:", guessed_db)


    df_proks = evaluate.generate_df_zymo()
    # Generate set of taxa names that belongs to the Zymo at the given level:
    set_proks = set(map(lambda a_str: a_str[3:], df_proks[taxo_cutoff]))
    # If 'species', str for genus need to be concatenated to the one for species
    if taxo_cutoff == 'species': 
        list_s = map(lambda a_str: a_str[3:], df_proks[taxo_cutoff])
        list_g = map(lambda a_str: a_str[3:], df_proks['genus'])
        set_proks = set([tupl[0] + " " + tupl[1] 
                         for tupl in zip(list_g, list_s)])


    # Guess the "mode":
    # Mode for handling of reads-clustering results
    CLUST_MODE = "_to" not in infile_base 
    IS_SAM_FILE = ext_infile == ".sam"

    if not CLUST_MODE and not IS_SAM_FILE:
        print("CENTRIFUGE MODE !\n")
        dict_gethered = {}

        with open(to_infile, 'r') as in_tab_file:
            in_tab_file.readline() # Skip header
            set_ignored_hits = set()
            set_all_readIDs = set()

            for line in in_tab_file:
                idx_first_tab = line.find('\t')
                readID = line[0:idx_first_tab] 
                rest = line[(idx_first_tab+1): ].rstrip('\n')
                set_all_readIDs.add(readID)

                # /!\ We ignore hits with taxid=0 and ref_name='no rank':
                if rest.split('\t')[0:2] == ["no rank", "0"]:
                    set_ignored_hits.add(readID + '\t' + rest)
                else:
                    if readID not in dict_gethered.keys():
                        dict_gethered[readID] = [rest]
                    else:
                        dict_gethered[readID].append(rest)

        # Be sure to don't exclude any read:
        assert(len(set_all_readIDs) == len(dict_gethered.keys()))
        print("NB OF IGNORED HITS:", len(set_ignored_hits))

        readIDs, list_val, list_indexes = dict_gethered.keys(), [], []
        dict_problems = {'2071623':'37482', '585494':'573', '595593':'656366'}

        for readID in readIDs:
            list_indexes.append(readID)
            if len(dict_gethered[readID]) > 1:
                list_taxids = list(map(lambda a_str: a_str.split('\t')[1], 
                                   dict_gethered[readID]))
                if '0' in list_taxids:
                    print(readID);sys.exit()
                nb_trashes = sum(map(pll.is_trash, list_taxids))
                # nb_trashes = 0
                lineage = 's'.join(list_taxids)
                # In the case of multiple hits, all hits have the same score, 
                # 2ndBestScore etc, except for the hitLength, that can differ
                (_, _, score, _, _, queryLength, 
                 numMatches) = dict_gethered[readID][0].split('\t')
                list_hitLength = list(map(lambda a_str: a_str.split('\t')[4], 
                                      dict_gethered[readID]))
                if len(set(list_taxids)) == 1:
                    type_align = 'second_uniq'
                else:
                    type_align = 'second_plural'
                hitLength = 'm'.join(list_hitLength)
                secondBestScore = pd.np.nan

            else:
                tupl = dict_gethered[readID][0].split('\t')
                ( ref_name, taxID, score, secondBestScore,
                  pre_hitLength, queryLength, numMatches ) = tupl

                if taxID in dict_problems.keys():
                    taxID = dict_problems[taxID]

                if ref_name == 'unclassified':
                    type_align = 'unmapped'
                    (lineage, score, secondBestScore, hitLength, queryLength,
                     numMatches, nb_trashes) = [pd.np.nan] * 7
                else:
                    type_align = 'normal'    
                    lineage = ';' + taxID + ';'
                    hitLength = ';' + pre_hitLength + ';'
                    nb_trashes = int(pll.is_trash(taxID))
                    # nb_trashes = 0

            list_val.append([type_align, lineage, nb_trashes, score, 
                             secondBestScore, hitLength, queryLength, 
                             numMatches])
        del readID

        tupl_columns = ('type_align', 'lineage', 'nb_trashes', 'score', 
                        'secondBestScore', 'hitLength', 'queryLength', 
                        'numMatches')
        my_csv = pd.DataFrame(data=None, columns=tupl_columns,
                              index=list_indexes)
        for i in range(len(tupl_columns)):
            # Cuz df.assign() takes 1 single karg:
            tmp_dict = {tupl_columns[i]:list(map(lambda sublist: sublist[i], 
                                                 list_val))}
            my_csv = my_csv.assign(**tmp_dict)
        del i, tmp_dict, list_val, list_indexes

        print(my_csv.type_align.value_counts())
        print()
        # # taxfoo = evaluate.taxfoo


    if not CLUST_MODE and IS_SAM_FILE:
        print("SAM MODE !")

        # Path to the needed 'seqid2taxid' file:
        to_seqid2taxid = to_dbs + "Centri_idxes/" + guessed_db + "/seqid2taxid"
        print()

        to_out_file = infile_base + ".csv"
        if not osp.isfile(to_out_file):
            # To make correspond operon number and taxid:
            dict_seqid2taxid = {}
            with open(to_seqid2taxid, 'r') as seqid2taxid_file:
                for line in seqid2taxid_file:
                    splitted_line = line.rstrip('\n').split('\t')
                    dict_seqid2taxid[splitted_line[0]] = int(splitted_line[1])
                del line

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
                    if query_name in dict_gethered.keys():
                    # if alignment.is_secondary or alignment.is_supplementary:
                        dict_gethered[query_name].append(dict_align)
                    else:
                        #assert(query_name not in dict_gethered.keys())
                        dict_gethered[query_name] = [dict_align]
                       
            input_samfile.close()
            del idx, alignment

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
                
                # Write mapped reads:
                for res_conversion in results:
                    to_write, list_header = str_from_res_conv(res_conversion)
                    if len(res_conversion) > 3:
                        # list_header = [a_key for a_key in res_conversion.keys() 
                        #                      if a_key != 'readID']
                        header_for_file = ','.join([''] + list_header)
                    out_file.write(to_write + '\n')
                del res_conversion
                # print(header_for_file);sys.exit()

                # Write unmapped reads:
                # if list_unmapped:
                if False:
                    for unmapped_read in list_unmapped:
                        out_file.write(unmapped_read + ',unmapped,no\n')
                    del unmapped_read

            # Write header:
            with open(to_out_file, 'r+') as out_file:
                content = out_file.read()
                out_file.seek(0)
                out_file.write(header_for_file + '\n' + content)
                # out_file.write(",type_align,lineage,nb_trashes,mapq,len_align," +
                #                "ratio_len,de\n")

            print("CSV CONVERSION TIME:", t.time() - TOTO)
            sys.exit()

                    
        else:
            print("FOUND CSV FILE:", to_out_file)
            print("Loading CSV file...")
            my_csv = pd.read_csv(to_out_file, header=0, index_col=0)
            print("CSV loaded !")

        # pd.Series(list_second_len, 
        #           name='Boxplot len_second (p1N300)').plot.box()#;plt.show()


    if not CLUST_MODE:
        print_distrib = False
        if print_distrib:
            test = my_csv[my_csv['lineage'].notnull()]['lineage']
            felix = pd.Series(map(lambda lin: len(lin.split('s')), test))
            bins_val = felix.max()
            felix.plot.hist(bins=bins_val, log=True)
            print("% DE + DE 25 SECOND:", sum(felix > 26)/len(felix)*100)
            # plt.plot([100]*bins_val)
            plt.show()
            # print(my_csv['type_align'].value_counts())
            # print(sum(my_csv['type_align']))
            sys.exit()

        TIME_CSV_TREATMENT = t.time()
        with_lineage = ((my_csv.type_align != 'unmapped') & 
                        (my_csv.type_align != 'only_suppl'))
                        # (my_csv.type_align != 'only_suppl'),
                        # (my_csv.type_align != 'second_plural')
        if not IS_SAM_FILE:
        # if False:
            cutoff_on_centriScore, cutoff_on_centriHitLength = 300, 50
            with_lineage = ((my_csv.type_align != 'unmapped') & 
                            (my_csv.type_align != 'only_suppl') &
                            (pd.to_numeric(my_csv['score']) > cutoff_on_centriScore))
            print("CUTOFF CENTRI SCORE:", cutoff_on_centriScore)
            print("NUMBER OF READS REMAINING:", sum(with_lineage), " | ",
                  "NB REMOVED:", len(with_lineage)-sum(with_lineage))

        my_csv_to_pll = my_csv[['lineage', 'type_align']][with_lineage]
        # MODE = 'LCA'
        MODE = 'MAJO'
        # if IS_SAM_FILE:
        #     MODE = 'MAJO'

        print()
        print("MODE FOR HANDLING OF THE MULTI-HITS =", MODE)

        partial_eval = partial(pll.eval_taxo, two_col_from_csv=my_csv_to_pll,
                                              set_levels_prok=set_proks,
                                              taxonomic_cutoff=taxo_cutoff,
                                              mode=MODE)
        print("PROCESSING CSV TO EVALUATE TAXO...")
        # Parallel version:
        eval_pool = mp.Pool(15)
        results_eval = eval_pool.map(partial_eval, my_csv_to_pll.index)
        eval_pool.close()
        # results_eval = list(map(partial_eval, my_csv_to_pll.index)) 
        del my_csv_to_pll
        print("PLL PROCESS FINISHED !")
        # sys.exit()
        
        # Compute counts:
        dict_count = {}
        counts_type_align = my_csv.type_align.value_counts()
        for type_aln, count in counts_type_align.items():
            dict_count[type_aln] = count
        del type_aln
        if 'second_uniq' not in counts_type_align.index:
            dict_count['second_uniq'] = 0
        if 'second_plural' not in counts_type_align.index:
            dict_count['second_plural'] = 0
        if 'unmapped' not in counts_type_align.index:
            dict_count['unmapped'] = 0
        dict_count['tot_second'] = (dict_count['second_plural'] + 
                                    dict_count['second_uniq'])
        dict_count["tot_reads"] = (dict_count['unmapped'] + 
                                   dict_count["tot_second"] +
                                   dict_count["normal"])        
        dict_stats['FN'] += dict_count['unmapped']

        # Add the results of the evaluation to the csv:
        dict_species2res = {} # To access evaluation results of a given species
        list_index_new_col, list_val_new_col = [], []
        for tupl_res in results_eval:
            readID, final_taxid, species, res_eval, remark_evaluation = tupl_res
            if remark_evaluation != 'no_majo_found':
                dict_species2res[species] = res_eval
                    
            list_index_new_col.append(readID)
            list_val_new_col.append((final_taxid, species, res_eval, 
                                     remark_evaluation))
        del tupl_res

        tmp_df = pd.DataFrame({'final_taxid':[tupl[0] 
                                          for tupl in list_val_new_col],
                               'species':[tupl[1] 
                                          for tupl in list_val_new_col],
                               'res_eval':[tupl[2] 
                                           for tupl in list_val_new_col],
                               'remark_eval':[tupl[3] 
                                              for tupl in list_val_new_col]},
                              index=list_index_new_col)
        my_csv = my_csv.assign(final_taxid=tmp_df.final_taxid,
                               species=tmp_df.species, 
                               res_eval=tmp_df.res_eval,
                               remark_eval=tmp_df.remark_eval)
        del tmp_df, list_val_new_col
        print("TIME FOR CSV PROCESSING:", t.time() - TIME_CSV_TREATMENT)
        print()

        # print(my_csv[pd.to_numeric(my_csv.score) <= 300]['res_eval'].value_counts())
        # sys.exit()

        
        # OTU mapping file writting:
        # (NaN values are automatically EXCLUDED during the 'groupby')
        grped_by_fin_taxid = my_csv.groupby(by=['final_taxid'])
        tool_used = 'centri'
        if IS_SAM_FILE:
            tool_used = 'last'
        sampl_prefix = (tool_used + guessed_db.capitalize() + 
                        MODE.lower().capitalize())

        with open(infile_base + '_' + MODE + '.map', 'w') as my_map:
            #my_map.write('#OTU ID\t' + 'SampleID' + '\n')
            for taxid_grp, grp in grped_by_fin_taxid:
                readIDs_to_write = map(lambda readID: sampl_prefix + '_' + 
                                                      str(readID), 
                                       grp.index)
                my_map.write(str(taxid_grp) + '\t' + 
                             '\t'.join(readIDs_to_write) + '\n')
                
                # taxo_to_write = ';'.join(['Other'] * 7)
                if taxid_grp != 'no_majo_found':
                    pass
                    # taxo_to_write = evaluate.taxo_from_taxid(taxid_grp)
                else:
                    print("NB OF 'NO_MAJO':", len(grp))
        print("--> Wrote OTUs mapping file for '{}' !".format(sampl_prefix))
        print()
        # sys.exit()


        # Count TP and FP for statistics:
        cutoff_nb_trashes = my_csv.nb_trashes.max()
        # cutoff_nb_trashes = 4
        print("CUTOFF ON THE NB OF TRASHES:", cutoff_nb_trashes)
        is_FP = ((my_csv.res_eval == 'FP') &
                 (my_csv.nb_trashes <= cutoff_nb_trashes))
        is_TP = ((my_csv.res_eval == 'TP') &
                 (my_csv.nb_trashes <= cutoff_nb_trashes))

        print("FP STATS:")
        # print(my_csv[is_FP][['remark_eval', 'nb_trashes']].groupby(['remark_eval', 'nb_trashes']).size())
        # print(my_csv[is_FP].remark_eval.value_counts().sort_index())
        print(my_csv[is_FP].nb_trashes.value_counts().sort_index())
        # print(my_csv[is_FP & (my_csv['remark_eval'] == 'no_majo_found')]['species'].value_counts())
        print()
        print("TP STATS:")
        print(my_csv[is_TP].species.value_counts())
        # print(my_csv[is_TP].type_align.value_counts().sort_index())
        # print(my_csv[is_TP][['remark_eval', 'nb_trashes']].groupby(['remark_eval', 'nb_trashes']).size())

        print()


        # test = my_csv[is_FP & (my_csv.final_taxid.notnull())]['final_taxid']
        # print(test.value_counts());sys.exit()
        # for readID, lin_val in test.items():
        #     # print([evaluate.taxfoo.get_taxid_name(int(taxid)) for taxid in lin_val.split('s')])
        #     my_lca = evaluate.taxfoo.find_lca(map(int, lin_val.split('s')))
        #     print(evaluate.taxfoo.get_taxid_rank(my_lca), evaluate.taxfoo.get_taxid_name(my_lca))
        

        # eval_taxfoo = evaluate.taxfoo
        # pd.DataFrame({'TP_nb_trashes':my_csv[is_TP]['nb_trashes'].value_counts(), 
        #               'FP_nb_trashes':my_csv[is_FP]['nb_trashes'].value_counts()
        #               }).plot.bar(title="Distrib nb_trashes between FP and TP",
        #                           color=('red', 'green'))
        # plt.show()

        # if not IS_SAM_FILE:
            # tupl_columns = ('type_align', 'lineage', 'nb_trashes', 'score', 
            #             'secondBestScore', 'hitLength', 'queryLength', 
            #             'numMatches')
            # my_series = pd.to_numeric(my_csv[is_TP]['score'])
            # my_series.plot.hist(bins=max(my_series), log=True)
            # plt.show()


        dict_stats['FP'] += sum(is_FP)
        dict_stats['FP'] += sum(my_csv['type_align'] == 'only_suppl')
        dict_stats['TP'] += sum(is_TP)
        # print(dict_stats)
        # print("TOT READS EVALUATED:", sum(dict_stats.values()))
        print()


        # Print general counting results:
        print(dict_count)
        tot_reads = dict_count['tot_reads']
        print("TOT READS:", tot_reads)
        print("% UNMAPPED:", round(dict_count["unmapped"]/tot_reads*100, 4))
        print()


        # AT THE TAXA LEVEL:
        print("RESULTS AT THE TAXA LEVEL:")
        print(dict_species2res)
        list_FP = [val for val in dict_species2res 
                   if dict_species2res[val] == 'FP']
        list_TP = [val for val in dict_species2res 
                   if dict_species2res[val] == 'TP']
        nb_pos_to_find = len(set_proks)
        recall_at_taxa_level = len(list_TP)/nb_pos_to_find
        print("TO FIND:", list_TP)
        print("NB_FP", len(list_FP), " | NB_TP", len(list_TP))
        dict_stats_sp_level = {'TP':len(list_TP),
                               'FP':len(list_FP)}
        precision_at_taxa_level = compute_metrics(dict_stats_sp_level, True)
        print("SENSITIVITY:", recall_at_taxa_level, " | ", "FNR:",
              1 - recall_at_taxa_level)
        calc_f1 = lambda tupl: round((2*tupl[0]*tupl[1])/(tupl[0]+tupl[1]), 4)
        print("F1-SCORE:", 
              calc_f1((precision_at_taxa_level, recall_at_taxa_level)))
        print()


        sum_stats, sum_counts = sum(dict_stats.values()),sum(dict_count.values())
        # assert(sum_stats == tot_reads - nb_evaluated_suppl)
        print("RESULTS AT THE READ LEVEL:")
        print(sorted(dict_stats.items()))
        dict_to_convert = dict_stats


    if CLUST_MODE:
        print('PAS LAAAAAA');sys.exit()
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
                res = evaluate.in_zymo(df_hits.loc[clust, "taxid"], taxo_cutoff, 
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
    print()
