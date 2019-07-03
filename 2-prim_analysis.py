#!/usr/bin/env python3


"""
Mapping statistics computation

Usage:
  2-prim_analysis.py (-i <inFile>) (-l <taxoCut>) [-w <writeMap>] [-o <outDir>]
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inFile=input_file     input file
  -l --taxoCut=taxo_cutoff   cutoff for the taxonomic level
  -w --writeMap=write_map    flag to either write (reads-vs-OTUs) mapping file [default: F] 
  -o --outDir=out_dir        /path/to/output_csv [default: ./]
"""


import sys, os
import os.path as osp
import time as t
import subprocess as sub
import pysam as pys
import multiprocessing as mp
import matplotlib.pyplot as plt
from itertools import islice
from functools import partial
from docopt import docopt
import src.check_args as check


def type_align_f(lin_val_from_df):
    """
    Tiny function used to determine the type of alignment, based on the 'lineage'
    value (string like: ';123;' or ';123s456;')
    """
    splitted_lineage = lin_val_from_df.strip(';').split('s')
    if len(splitted_lineage) == 1: # NORMAL
        return 'normal'
    else:
        set_lin = set(splitted_lineage)
        if len(set_lin) == 1: # SEVERAL IDENTIC TAXIDS
            return 'second_uniq'
        return 'second_plural'


def transform_ONT_CSV(to_SAM_csv):
    """
    Transform the CSV produced by Epi2me to make it usable by my scripts
    """
    base_filename = osp.basename(to_SAM_csv).lower()
    sep_used = ','

    if "epi2me" in base_filename:
        detected_tool = 'EPI2ME'
        col_to_take = ['read_id', 'taxid', 'exit_status']
        success = 'Classification successful'
    elif "wimp" in base_filename:
        detected_tool = 'WIMP'
        col_to_take = ['readid', 'taxID', 'exit_status']
        success = 'Classified'
    else: # Centrifuge-generated CSV
        detected_tool = 'Centrifuge'
        sep_used = '\t'
        col_to_take = ['readID', 'seqID', 'taxID', 'score', 'hitLength', 
                       'queryLength']

    print("Transforming CSV generated with {} ONT tool !".format(detected_tool))
    initial_csv = pd.read_csv(to_SAM_csv, header=0, sep=sep_used, 
                              usecols=col_to_take)

    # Homogeneization of columns names:
    col_conv = {'readid':'readID', 'read_id':'readID', 'taxid':'taxID'}
    initial_csv.columns = [col_conv.get(x, x) for x in initial_csv.columns]

    # Actually NO readID... :
    if all(initial_csv.readID.isnull()):
        print("NO readID actually found --> adding a pseudo-readID..")
        nb_rows = len(initial_csv.index)
        initial_csv.readID = [foo[0]+str(foo[1]) 
                                    for foo 
                                    in zip(['read_']*nb_rows, range(nb_rows))]

    if detected_tool == 'Centrifuge': # Filter on centriScore and hitLen here
        # First, we detect and remove odd entries:
        are_odd_entries = ((initial_csv.seqID == 'no rank') & 
                           (initial_csv.taxID == 0))
        print("NB OF IGNORED HITS ('no rank' + taxid=0):", sum(are_odd_entries))
        initial_csv = initial_csv[-are_odd_entries] # We filter them out

        # Then we treat 'unmapped':
        are_unass = initial_csv.seqID == 'unclassified'
        are_assigned = -are_unass

        # Then we apply a filter to limit the number of FP:
        # (In the case of multiple hits, all hits have the same score, 
        # 2ndBestScore etc, except for the hitLength, that can differ)
        cutoff_centriScore, cutoff_centriHitLen = 300, 50
        # cutoff_centriScore, cutoff_centriHitLen = 0, 0
        centri_filt = initial_csv.score > cutoff_centriScore
        centri_filt = (centri_filt & 
                       (initial_csv.hitLength > cutoff_centriHitLen))

        # KEY STEP:
        initial_csv = initial_csv[centri_filt]

        # def get_min_hitLength(hitLength_val):
        #     if type(hitLength_val) == str:
        #         list_val = hitLength_val.strip(';').split('m')
        #         if len(list_val) > 1:
        #             return min(map(float, list_val))
        #         return float(list_val[0])
        #     return hitLength_val

        # my_csv['min_hitLength'] = my_csv.hitLength.apply(get_min_hitLength)
        # centri_filter = my_csv.score > cutoff_on_centriScore
        # centri_filter = (centri_filter & 
        #                  (my_csv.min_hitLength > cutoff_on_centriHitLength))

        nb_nan = sum(are_unass)
        nb_pass_before_filter = sum(are_assigned)
        pass_after_filter = are_assigned & centri_filt
        nb_removed = nb_pass_before_filter-sum(pass_after_filter)
        assert(len(pass_after_filter) == sum(pass_after_filter)+nb_removed+nb_nan)

        nb_not_nan = len(pass_after_filter)-nb_nan
        if not cutoff_centriScore:
            print("   >> CENTRI NOT FILTERED")
        else:
            print("   >> CUTOFF CENTRI SCORE:", cutoff_centriScore, 
                  " | ", "CUTOFF CENTRI HitLength:", 
                  cutoff_centriHitLen)
        print("NB OF UNCLASSIF:", nb_nan, 
              ' | NB_NOT_NaNs:', nb_not_nan)
        print("NB REMOVED (NaNs excepted):", nb_removed, 
              "~ {}% removed".format(int(nb_removed/nb_not_nan*100)))

        # Index has changed, so we need to redefine these variables:
        are_unass = initial_csv.seqID == 'unclassified'
        are_assigned = -are_unass
    
    else:
        # Treat 'unmapped' here too:
        # /!\ This func assert that every status different from 
        # 'Classif successful' can be considered as 'unmapped'
        are_assigned = initial_csv.exit_status == success
        are_unass = -are_assigned

    nb_unass = sum(are_unass)
    if detected_tool != 'Centrifuge':
        print("NB UNASS:", nb_unass)

    print("Formatting properly (for further taxo eval)..")
    # Only assigned reads are taken for the next part:
    initial_assigned = initial_csv[are_assigned]
    
    # if len(initial_csv[name_col_readID]) == len(initial_csv[name_col_readID].unique()):

    # Solve problematic taxids (mostly with 'p_compressed'):
    dict_problems = {'2071623':'37482', '585494':'573', '595593':'656366',
                     '697046':'645', '1870930':'1812935', 
                     '2036817':'106590', '134962':'53431', 
                     '1662457':'553814', '1662456':'553814',
                     '1662458':'553814', '1834200':'1796646',
                     '1902251':'745377', '913102':'59277'}
    initial_csv.loc[are_assigned , 'taxID'] = initial_assigned.taxID.apply(
                                lambda val: dict_problems.get(str(val), val))

    # Tiny functions used to obtain needed values (str lineage and nb trashes):
    lineage_f = lambda a_arr: (';' + 's'.join(
                                        str(taxid) for taxid in a_arr) + ';')
    nb_trash_f = lambda a_arr: sum(map(pll.is_trash, a_arr))

    # Group by read and determine needed values (use of Numpy to speed it up):
    keys, values = initial_csv[['readID', 'taxID']][are_assigned].sort_values('readID').values.T
    ukeys, index = pd.np.unique(keys, return_index=True)
    arrays = pd.np.split(values, index[1:])
    grped_csv = pd.DataFrame({'readID':ukeys, 
                              'lineage':[lineage_f(a) for a in arrays],
                              'nb_trashes':[nb_trash_f(a) for a in arrays]})

    # Add 'type_align' column and 'unmapped' reads:
    grped_csv['type_align'] = grped_csv.lineage.apply(type_align_f)
    tmp_df = pd.DataFrame({'type_align':['unmapped'] * nb_unass}, 
                          index=initial_csv[are_unass].readID)
    return pd.concat([grped_csv.set_index('readID'), tmp_df], sort=True)


def ali_to_dict(align_obj):
    """
    Takes an alignment instance and return a dict containing the needed
    attributes (only)

    With infer_read_length() method, hard-clipped bases are included in the 
    counting
    infer_query_length() method does NOT include them
    """
    if align_obj.is_unmapped:
        return 'unmapped'

    ratio_len = align_obj.infer_query_length()/align_obj.infer_read_length()
    assert(ratio_len <= 1)

    to_return = {"mapq":align_obj.mapping_quality,
                 "ref_name": align_obj.reference_name,
                 "len_align":align_obj.infer_query_length(),
                 "ratio_len":ratio_len,
                 "is_suppl":align_obj.is_supplementary,
                 "has_SA": align_obj.has_tag("SA"),
                 "is_second":align_obj.is_secondary}
                 
    if align_obj.has_tag('de'):
        to_return['de'] = round(align_obj.get_tag("de"), 4)
    if align_obj.has_tag('AS'):
        to_return['AS'] = align_obj.get_tag("AS")

    return to_return


def str_from_res_conv(dict_res_conv):
    """
    Generate a string (ready to write) from a SAM to CSV conversion
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

    return (",".join(to_return + [lineage_to_write] + 
                     [str(x[1]) for x in rest]), 
            to_extract[1:] + list(map(lambda x: x[0], rest)))


def get_ancester_name(arg_taxid_ancest):
    """
    Return the sp_name of the ancester given its taxid (different of what I 
    called 'final_taxid')
    If a taxid could not be found for the ancester, 'notDeterminable' is 
    returned
    """
    if arg_taxid_ancest == 'no_majo_found':
        return 'notDeterminable'
        
    if not pd.np.isnan(arg_taxid_ancest):
        return taxfoo.get_taxid_name(arg_taxid_ancest)
    return 'notDeterminable'


def discriminate_miss(arg_taxid, wanted_taxo, df_proks_arg):
    """
    Given the taxid of a miss, discriminate between FFP (= in Zymo but taxo too 
    high) and TFP (= 'final_taxid' misassigned or taxo lvl above 'Bacteria')
    """
    lineage = taxfoo.get_dict_lineage_as_taxids(arg_taxid, 
                                                want_taxonomy=wanted_taxo)
    for taxo_lvl in wanted_taxo:
        if taxo_lvl in lineage.keys():
            current_ancester = lineage[taxo_lvl]
            current_set_proks = set(map(lambda a_str: a_str[3:], 
                                        df_proks_arg[taxo_lvl]))
            res_eval = evaluate.in_zymo(current_ancester, 
                                        current_set_proks, taxo_lvl)[-1]
            if res_eval == 'true_pos':
                return pd.Series(['FFP', current_ancester,
                                  taxfoo.get_taxid_name(current_ancester)])

    del taxo_lvl

    return pd.Series(['TFP', arg_taxid, taxfoo.get_taxid_name(int(arg_taxid))])


def dict_stats_to_vectors(dict_res):
    """
    Generate 2 vectors (for predicted and true values), based on the values of
    'well', 'miss' etc contained in the dict_res given as arg
    """
    vec_pred, vec_true = [], []

    for res in dict_res:
        if res == 'well':
            vec_true.extend([True]*dict_res[res]) # Positive in reality
            vec_pred.extend([True]*dict_res[res]) # Well predicted positive
        elif res == 'unass':
            vec_true.extend([True]*dict_res[res]) # Positive in reality
            vec_pred.extend([False]*dict_res[res]) # But predicted negative
        else: # if res == 'miss':
            vec_true.extend([False]*dict_res[res]) # Negative in reality
            vec_pred.extend([True]*dict_res[res]) # But predicted positive

    return (vec_true, vec_pred)


def compute_metrics(dict_stats_to_convert, at_taxa_level):
    """
    Compute different metrics with sklearn, from a dict counting miss, well and
    unass
    """
    import sklearn.metrics as skm
    y_true, y_pred = dict_stats_to_vectors(dict_stats_to_convert)

    precision = round(skm.precision_score(y_true, y_pred), 4) # PPV = well/(well+miss) (positive predictive value or precision)
    print("PRECISION:", precision, " | ", "FDR:", round(1 - precision, 4))

    if not at_taxa_level:
        # sensitivity = round(skm.recall_score(y_true, y_pred), 4)
        tp, fp = dict_stats_to_convert['well'], dict_stats_to_convert['miss'] # CORRECTION !
        sensitivity = round(tp/(tp + fp + dict_stats_to_convert['unass']), 4)
    
        print("SENSITIVITY:", sensitivity, " | ", "FNR:", 
              round(1 - sensitivity, 4)) # TPR = well/(well+unass) (or sensitivity, hit rate, recall) 
    else:
        return precision



# MAIN:
if __name__ == "__main__":
    # START OF ARGUMENTS CHECKING:
    ARGS = docopt(__doc__, version='0.1')
    to_infile, infile_base, _, _, ext_infile = check.infile(ARGS["--inFile"], 
                                                ['csv', 'tsv', 'sam'])
    write_map = check.bool_type(ARGS["--writeMap"])
    outDir = ARGS["--outDir"]
    assert(osp.isdir(outDir))

    # Common variables:
    NB_THREADS = 10
    dict_stats = {'unass':0, 'well':0, 'miss':0}

    print()
    print("Loading taxonomic Python module...")
    import src.parallelized as pll
    evaluate = pll.eval
    pd = evaluate.pd
    taxfoo = evaluate.taxfoo
    print("Taxonomic Python module loaded !\n")
    taxo_cutoff = check.acceptable_str(ARGS["--taxoCut"], evaluate.want_taxo)

    print("PROCESSING FILE:", infile_base)
    print("OutDir:", outDir)
    print("Script PARALLELIZED on {} CPUs".format(NB_THREADS))

    # Detect the tool used:
    IS_SAM_FILE = ext_infile == ".sam"
    lower_infile_base = infile_base.lower()

    if IS_SAM_FILE:
        tool_used = 'minimap2'
    else:
        if 'centri' in lower_infile_base:
            tool_used = 'centri'    
        elif 'epi2me' in lower_infile_base:
            tool_used = 'epi2me'
        elif 'wimp' in lower_infile_base:
            tool_used = 'wimp'
        else:
            print("ERROR: With a CSV (or TSV) file, the tool used to " + 
                  "generate it must be specified within the filename") 
            print("(ideally like that: 'epi2me_your_fav_file.csv')")
            print() ; sys.exit()
    print("DETECTED TOOL:", tool_used)

    # Guess the database used:
    is_ONT_tool = tool_used in ('epi2me', 'wimp')
    pos_name = infile_base.find("_to")

    if is_ONT_tool: # NO '_to' needed in this case
        guessed_db = 'ONTcomm'
    elif pos_name == -1: # '_to' not found --> PROBLEM
        print("Database name NOT specified within filename " + 
              "(e.g. 'your_name_toRrn.sam')")
        print() ; sys.exit(2)
    else:
        # /!\ This step is problematic if DB_name contains an underscore:
        guessed_db = infile_base[(pos_name+3):].split('_')[0]

    print("DB GUESSED:", guessed_db) ; print()


    # Get 'seqid2taxid' path from 'pipeline.conf', in case of Minimap:
    if IS_SAM_FILE:
        to_conf_file = 'pipeline.conf'
        conf_csv = pd.read_csv(to_conf_file, sep=';', comment='#')
        # Lowercase conversion, to make it more flexible:
        conf_csv.tool.str.lower() ; conf_csv.type_param.str.lower()
        params_csv = conf_csv[conf_csv.tool == 'minimap2'].drop(['tool'], axis='columns')
        cond_to_seqid2taxid = ((params_csv.type_param == 'to_seqid2taxid') & 
                               (params_csv.supplField_1.str.lower() == guessed_db.lower()))
        list_to_seqid2taxid = params_csv[cond_to_seqid2taxid].supplField_2.values
        to_seqid2taxid = check.list_config(list_to_seqid2taxid, 
                                           'to_seqid2taxid')
        print("Path to the 'seqid2taxid' file deducted from", to_conf_file)


    # END of args checking --> The core program can START:

    # Generate set of taxa names that belongs to the Zymo at the given level:
    df_proks = evaluate.generate_df_zymo()
    set_proks = set(map(lambda a_str: a_str[3:], df_proks[taxo_cutoff]))
    # If 'species', str for genus need to be concatenated to the one for species
    if taxo_cutoff == 'species': 
        list_s = map(lambda a_str: a_str[3:], df_proks[taxo_cutoff])
        list_g = map(lambda a_str: a_str[3:], df_proks['genus'])
        set_proks = set([tupl[0] + " " + tupl[1] 
                         for tupl in zip(list_g, list_s)])

    
    # We start the parsing, according to the type of input file:
    if not IS_SAM_FILE:
        # If CSV file as input (Centrifuge or ONT tools WIMP/EPI2ME), we just 
        # need to reformat this CSV to make it usable:
        formatting_start = t.time()
        my_csv = transform_ONT_CSV(to_infile)
        print(my_csv.type_align.value_counts())
        print()
        print("ELAPSED TIME FOR FORMATTING:", t.time()-formatting_start)
        # sys.exit()


    else: # IS_SAM_FILE == TRUE:
        print("SAM MODE !")
        print()

        to_out_file = osp.join(outDir, infile_base + ".csv")
        if not osp.isfile(to_out_file):
            # To make correspond sequence id (fasta header) and taxid:
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

                else:
                    dict_align = ali_to_dict(alignment)
                    if query_name in dict_gethered.keys():
                        dict_gethered[query_name].append(dict_align)
                    else:
                        dict_gethered[query_name] = [dict_align]
                       
            input_samfile.close()
            del idx, alignment

            print("SAM PARSING TIME:", str(t.time() - START_SAM_PARSING))
            print()
            assert(len(list_unmapped) == len(set(list_unmapped)))


            # Extract relevant SAM info towards a CSV file:
            TOTO = t.time()
            print("Converting SAM into CSV...")
            conv_partial = partial(pll.SAM_to_CSV, 
                                   conv_seqid2taxid=dict_seqid2taxid)


            # VERSION 2 (use much less RAM but a bit slower):
            dict_light = (tupl for tupl in dict_gethered.items())
            nb_reads_to_process = len(dict_gethered.keys())
            proc_chunks, chunksize = [], nb_reads_to_process//NB_THREADS

            while True:
                my_slice = list(islice(dict_light, chunksize))
                if my_slice:
                    proc_chunks.append(my_slice)
                else:
                    break
            assert(sum(map(len, proc_chunks)) == nb_reads_to_process)

            def process_chunk(chunk, partial_fct):
                return list(map(partial_fct, chunk))

            results = []
            with mp.Pool(NB_THREADS) as my_pool:
                proc_results = [my_pool.apply_async(process_chunk, 
                                                    args=(chunk, conv_partial))
                                for chunk in proc_chunks]
                # Magic line that if basically a one-liner of 'extend' method:
                results = [r for r_ext in proc_results for r in r_ext.get()]

            print("PLL PROCESS FINISHED !")
            dict_gethered.clear()


            # Write outfile:
            with open(to_out_file, 'w') as out_file:
                # Write mapped reads:
                for dict_res_conv in results:
                    to_write, list_header = str_from_res_conv(dict_res_conv)
                    if len(dict_res_conv) > 3:
                        header_for_file = ','.join([''] + list_header)
                    out_file.write(to_write + '\n')
                del dict_res_conv

                # Write unmapped reads:
                if list_unmapped:
                    for unmapped_read in list_unmapped:
                        out_file.write(unmapped_read + ',unmapped,no\n')
                    del unmapped_read

            # Add header to the outfile:
            with open(to_out_file, 'r+') as out_file:
                content = out_file.read()
                out_file.seek(0)
                out_file.write(header_for_file + '\n' + content)

            print("CSV CONVERSION TIME:", t.time() - TOTO)
            # sys.exit()
            my_csv = pd.read_csv(to_out_file, header=0, index_col=0)

                    
        else:
            print("FOUND CSV FILE:", to_out_file)
            print("Loading CSV file...")
            my_csv = pd.read_csv(to_out_file, header=0, index_col=0)
            print("CSV loaded !")



    if True:
        # To print the distribution of the number of equivalent hits by reads:
        print_distrib = False
        if print_distrib:
            test = my_csv.lineage.dropna()
            felix = pd.Series(map(lambda lin: len(lin.split('s')), test))
            bins_val = felix.max()
            felix.plot.hist(bins=bins_val, log=True)
            print("% DE + DE 25 SECOND:", sum(felix > 26)/len(felix)*100)
            plt.show()
            sys.exit()


        TIME_CSV_TREATMENT = t.time()
        with_lineage = ((my_csv.type_align != 'unmapped') & 
                        (my_csv.type_align != 'only_suppl'))

        # Filter only reads that are taxonomically evaluable (assigned):
        my_csv_to_pll = my_csv[['lineage', 'type_align']][with_lineage].reset_index()
        nb_reads_to_process = len(my_csv_to_pll.index)


        # Compute counts:
        dict_count = {'second_uniq':0, 'second_plural':0, 'unmapped':0,
                      'only_suppl':0}
        counts_type_align = my_csv.type_align.value_counts()
        for type_aln, count in counts_type_align.items():
            dict_count[type_aln] = count
        del type_aln, count
        dict_count['tot_second'] = (dict_count['second_plural'] + 
                                    dict_count['second_uniq'])
        dict_count["tot_reads"] = sum(counts_type_align)
        dict_count["tot_mapped"] = nb_reads_to_process 
        dict_stats['unass'] += dict_count['unmapped']

        # Print general counting results:
        print()
        print(dict_count)
        tot_reads = dict_count['tot_reads']
        print("TOT READS:", tot_reads)
        print("% UNMAPPED:", round(dict_count["unmapped"]/tot_reads*100, 2))
        

        MODE = 'LCA'
        MODE = 'MAJO_OLD'
        # MODE = 'TOP_ONE'
        MODE = 'MINOR_RM_LCA'

        print("\nMODE FOR HANDLING MULTI-HITS:", MODE)
        print()
        print("PROCESSING {} reads from CSV to EVALUATE TAXO..".format(
                                                        nb_reads_to_process))
        print("(Tot nb of entries: {})".format(
                                sum(map(lambda lin: len(lin.split('s')), 
                                        my_csv_to_pll.lineage.values))))

        partial_eval = partial(pll.eval_taxo, set_levels_prok=set_proks,
                                              taxonomic_cutoff=taxo_cutoff,
                                              mode=MODE)

        # Parallel version:
        def process_chunk(chunk_df, partial_func):
            return chunk_df.apply(partial_func, axis='columns')


        proc_chunks, chunksize = [], nb_reads_to_process//NB_THREADS
        if nb_reads_to_process < NB_THREADS: # Case when the nb of reads is too small
            chunksize = 1

        for i_proc in range(NB_THREADS):
            chunkstart = i_proc * chunksize
            # make sure to include the division remainder for the last process
            chunkend = (i_proc + 1) * chunksize if i_proc < NB_THREADS - 1 else None
            proc_chunks.append(my_csv_to_pll.iloc[chunkstart : chunkend])
        assert(sum(map(len, proc_chunks)) == nb_reads_to_process)

        with mp.Pool(NB_THREADS) as eval_pool:
            proc_results = [eval_pool.apply_async(process_chunk, 
                                                  args=(chunk, partial_eval))
                            for chunk in proc_chunks]
            result_chunks = [r.get() for r in proc_results]
        print("PLL PROCESS FINISHED !")

        del my_csv_to_pll

        print('Finalizing evaluation..')
        # Add res to the main CSV:
        my_res = pd.concat(result_chunks, sort=True)
        col_to_add = ['taxid_ancester', 'final_taxid', 'res_eval', 
                      'remark_eval']
        my_csv = pd.concat([my_csv, my_res.set_index('index')[col_to_add]], 
                           axis='columns', sort=True, copy=False)
        list_col_my_res = my_res.set_index('index').columns.values.tolist()
        list_col_my_csv = my_csv.columns.values.tolist()
        assert(sorted(my_csv.columns) == sorted(set(list_col_my_res+list_col_my_csv)))


        # CORRECTION of B. intestinalis:
        if guessed_db.lower() == 'rrn':
            print("\n  >>> CORRECTION OF INTESTINALIS FOR RRN!")
            are_intestinalis = my_csv.final_taxid == 1963032
            print("(NB of corrected reads: {})".format(sum(are_intestinalis)))
            print()
            my_csv.loc[are_intestinalis, 'res_eval'] = 'well'
            my_csv.loc[are_intestinalis, 'final_taxid'] = 1423
            ancest_subtilis = taxfoo.get_dict_lineage_as_taxids(1423)[taxo_cutoff]
            my_csv.loc[are_intestinalis, 'taxid_ancester'] = ancest_subtilis


        my_csv['species'] = my_csv.taxid_ancester.apply(get_ancester_name)


        # To access evaluation results of a given species:
        tmp_df = my_csv[['species', 'res_eval']].drop_duplicates()
        dict_species2res = tmp_df.set_index('species')['res_eval'].to_dict()
        del tmp_df

        print("TIME FOR CSV PROCESSING:", 
              round(t.time() - TIME_CSV_TREATMENT, 4))
        print("NB of 'no_majo_found':", 
              sum(my_csv.final_taxid.astype(str) == 'no_majo_found'))
        print()

        
        # Mapping file writting (reads vs OTUs):
        if write_map:
            # (NaN values are automatically EXCLUDED during the 'groupby')
            grped_by_fin_taxid = my_csv.groupby(by=['final_taxid'])

            if IS_SAM_FILE:
                run_name = '.'.join(infile_base[0:pos_name].split('_'))
            else:
                index_underscore = infile_base.find('_') + 1
                run_name = '.'.join(infile_base[index_underscore:pos_name].split('_'))

            sampl_prefix = (run_name + '.' + tool_used + 
                            guessed_db[0].upper() + guessed_db[1:] + 
                            MODE.split('_')[0].capitalize())

            to_out_map = osp.join(outDir, infile_base + '_' + MODE + '.map')
            with open(to_out_map, 'w') as my_map:
                for taxid_grp, grp in grped_by_fin_taxid:
                    readIDs_to_write = map(lambda readID: sampl_prefix + '_' + 
                                                          str(readID), 
                                           grp.index)
                                        
                    if taxid_grp != 'no_majo_found':
                        # 'int' casting needed for the taxid
                        taxid_to_write = int(taxid_grp)
                    else:
                        taxid_to_write = taxid_grp
                        # Like this, 'no_majo_found' are excluded ? 
                        # continue 

                    my_map.write(str(taxid_to_write) + '\t' + 
                                 '\t'.join(readIDs_to_write) + '\n')

                unmapped_to_write = map(lambda readID: sampl_prefix + '_' + 
                                                          str(readID), 
                                        my_csv[my_csv.type_align == 'unmapped'].index)     
                my_map.write('unmapped\t' +  '\t'.join(unmapped_to_write) + 
                             '\n')
            print("  >> Wrote OTUs mapping file ({})".format(my_map.name))
            print("  >> With '{}' as sample name".format(sampl_prefix))
            print()


        # Count well and miss for statistics:
        cutoff_nb_trashes = my_csv.nb_trashes.max()
        # cutoff_nb_trashes = 4
        print("CUTOFF ON THE NB OF TRASHES:", cutoff_nb_trashes)
        is_miss = ((my_csv.res_eval == 'miss') &
                 (my_csv.nb_trashes <= cutoff_nb_trashes))
        is_well = ((my_csv.res_eval == 'well') &
                 (my_csv.nb_trashes <= cutoff_nb_trashes))

        print("WELL STATS:")
        print(my_csv[is_well].species.value_counts())
        print()


        # EXTRACTION of aligned sequeces within reference DB (use for custom_NanoSim):
        extract_ref = False
        if extract_ref:
            print("EXTRACTIING..")

            with pys.AlignmentFile(to_infile, "r") as input_samfile:
                dict_felix = {}
                for idx, alignment in enumerate(input_samfile.fetch(until_eof=True)):
                    if alignment.is_unmapped:
                        continue
                    elif alignment.has_tag("SA"):
                        continue
                    elif alignment.is_secondary:
                        continue
                    else:
                        query_name = alignment.query_name
                        # Take only well reads (miss suposed noisy):
                        if my_csv.loc[query_name, 'res_eval'] == 'well':
                            assert(alignment.query_name not in dict_felix.keys())
                            tmp_list = [alignment.reference_name, 
                                        alignment.reference_start, 
                                        alignment.reference_end]
                            dict_felix[alignment.query_name] = tmp_list
                del idx, alignment

            print("{} extractable sequences (miss only)".format(len(dict_felix)))
            test_df = pd.DataFrame.from_dict(dict_felix, orient='index', 
                                             columns=['ref_name', 'start_pos', 
                                                      'end_pos'])
            dict_felix.clear()
            print(test_df.head())
            test_df.to_csv('extractable_' + infile_base + '.csv', header=True)
            print("Done !")
            sys.exit()


        print("MISS STATS:")
        # print(my_csv[is_miss].remark_eval.value_counts())
        print(my_csv[is_miss & (my_csv.remark_eval == 'second_uniq;misassigned')].final_taxid.apply(taxfoo.get_taxid_name).value_counts())
        print()

        # Add numbers of well and miss to the dict of stats:
        dict_stats['miss'] += sum(is_miss)
        dict_stats['miss'] += sum(my_csv['type_align'] == 'only_suppl')
        dict_stats['well'] += sum(is_well)


        # AT THE TAXA LEVEL:
        print("RESULTS AT THE TAXA LEVEL:")
        list_miss = [val for val in dict_species2res 
                                 if dict_species2res[val] == 'miss']
        list_well = [val for val in dict_species2res 
                                 if dict_species2res[val] == 'well']
        nb_positive_to_find = len(set_proks)
        nb_well_taxa, nb_miss_taxa = len(list_well), len(list_miss)

        sens_taxa_lvl = len(list_well)/nb_positive_to_find
        print("TO FIND:", sorted(set_proks))
        print("NB_MISS", nb_miss_taxa, " | NB_WELL", nb_well_taxa)
        dict_stats_sp_level = {'well':nb_well_taxa,
                               'miss':nb_miss_taxa}
        prec_taxa_lvl = round(nb_well_taxa / (nb_well_taxa + nb_miss_taxa), 4)
        FDR_taxa_lvl = round(1 - prec_taxa_lvl, 4)
        print("PRECISION:  {} | FDR: {}".format(prec_taxa_lvl, FDR_taxa_lvl))
        print("SENSITIVITY:", sens_taxa_lvl, " | ", "FNR:", 1-sens_taxa_lvl)
        calc_f1 = lambda tupl: round((2*tupl[0]*tupl[1])/(tupl[0]+tupl[1]), 4)
        print("F1-SCORE:", calc_f1((prec_taxa_lvl, sens_taxa_lvl)))
        print()


        # Make difference between TFP and FFP:
        if False:
            print("DISCRIMINATION between TFP and FFP (at read lvl):")
            possible_notInKeys = my_csv.remark_eval.apply(
                                        lambda val: "notInKeys" in str(val))
            miss_notInKey = my_csv[is_miss & possible_notInKeys]

            if miss_notInKey.empty:
                print("NOT POSSIBLE ! (no miss and/or no 'notInKeys')")
            else:
                counts_miss_NotInKey = miss_notInKey.final_taxid.value_counts()
                counts_miss_NotInKey.name = 'counts'

                # We remove 'superkingdom' and 'species' lvls
                taxo_not_bact = evaluate.want_taxo[1:][::-1][1:] 
                df_miss = pd.DataFrame(counts_miss_NotInKey).reset_index()
                tmp_df = df_miss['index'].apply(discriminate_miss, 
                                              args=(taxo_not_bact, df_proks)) 
                df_miss = df_miss.assign(status=tmp_df[0], ancester_taxid=tmp_df[1],
                                     ancester_name=tmp_df[2])
                del tmp_df

                tot_Fmiss = sum(df_miss[df_miss.status == 'Fmiss'].counts)

                print('TOT NB OF Fmiss: {} | OVER {} miss_notInKey'.format(tot_Fmiss,
                                                                       len(miss_notInKey)))
                tot_miss = dict_stats['miss']
                print("Total of {} miss --> {} + {} othermiss".format(tot_miss, tot_Fmiss, 
                                                                  tot_miss-tot_Fmiss))
            print()


        sum_stats, sum_counts = sum(dict_stats.values()),sum(dict_count.values())
        print("RESULTS AT THE READ LEVEL:")
        print(sorted(dict_stats.items()))
        dict_to_convert = dict_stats

       
    # COMPUTE METRICS:
    compute_metrics(dict_to_convert, False)
    print()
