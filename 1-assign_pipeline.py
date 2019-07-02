#!/usr/bin/env python3

"""
Long-reads (LR) classification pipeline

Usage:
  1-assign_pipeline.py (-i <inputFqFile>) [-d <db>] [-T <trim>] [-C <chim>] [-D <deter>] [-t <threads>] [-o <outDir>]
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inputFqFile=input_fq  input fastq file
  -d --db=database           name of the database to use [default: Zymo]
  -T --trim=trimming         use Porechop to trim reads [default: no]
  -C --chim=chimera          chimera detection using yacrd [default: no]
  -D --deter=determination   taxonomic determination step [default: no]
  -t --threads=nb_threads    number of threads for parallelisation [default: 10]
  -o --outDir=out_dir        output directory for taxo file (SAM or CSV) [default: ./]
  
Arguments:
    trimming: 'porechop' or 'no' (default)
    chimera: 'yacrd' or 'no' (default)
    determination: 'minimap2', 'centri', or 'no' (default)
    database: 'SILVA', 'rrn', 'p_compressed', 'zymo' or 'newLot_Zymo'
"""

import os
import sys
import subprocess as sub
import time as t
import pandas as pd
import os.path as osp
from docopt import docopt
import src.check_args as check
        
        

# MAIN:
if __name__ == "__main__":
    START_TIME = t.time()
    
    # Get arguments:
    ARGS = docopt(__doc__, version='0.1')
    tuple_check_fq = check.infile(ARGS["--inputFqFile"], 
                                  ('fq', 'fastq', 'fa', 'fasta', 'mfasta'))
    in_fq_path, in_fq_base, _, tail_input_fq, in_fq_ext = tuple_check_fq
    DB_NAME = check.acceptable_str(ARGS["--db"], 
        ['rrn', 'zymo', 'newlot_zymo', 'silva', 'p_compressed', 'ncbi_16s'])
    TRIM = check.acceptable_str(ARGS["--trim"], ["porechop", "no"])
    CHIM = check.acceptable_str(ARGS["--chim"], ["yacrd", "no"])
    DETER = check.acceptable_str(ARGS["--deter"], 
                                 ["minimap2", "centri", "no"])
    nb_threads = check.input_nb(ARGS["--threads"], "'-t number of threads'")
    outDir = ARGS["--outDir"]
    if not osp.isdir(outDir):
        print("ERROR: outDir does NOT exit !") ; sys.exit()

    print()
    print("Tool: {} | Against: {}".format(DETER, DB_NAME))
    print("Output directory:", outDir)


    # Extracting parameters from 'pipeline.conf' file:
    conf_csv = pd.read_csv('pipeline.conf', sep=';', comment='#')
    # Lowercase conversion, to make it more flexible:
    conf_csv.tool.str.lower() ; conf_csv.type_param.str.lower()
    params_csv = conf_csv[conf_csv.tool == DETER].drop(['tool'], axis='columns')

    cond_to_tool = params_csv.type_param == 'to_exe'
    list_to_tool = params_csv[cond_to_tool].supplField_1.values
    if len(list_to_tool) > 1:
        print("ERROR: Several possible params for 'to_exe'") ; sys.exit()

    cond_to_ref_db = ((params_csv.type_param == 'to_ref') & 
                      (params_csv.supplField_1.str.lower() == DB_NAME))
    list_to_ref_db = params_csv[cond_to_ref_db].supplField_2.values
    to_ref_db = check.list_config(list_to_ref_db, 'to_ref')
    if DETER != 'centri':
        assert(osp.isfile(to_ref_db))
    

    print()
    print("Deducted from conf file:")
    print("To reference database:", to_ref_db)
    
                
    #Common variables/params:
    # To databases directory:
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    path_apps = "/home/sheldon/Applications/"
    path_proj = "/projets/metage_ONT_2019/"
    
    
    # Adaptators trimming using Porechop:
    path_to_porechop = "home/sheldon/Applications/Porechop-0.2.4/porechop-runner.py"
    args_porechop = " --discard_middle"
    dirOut_porechop = path_proj + "1_trim/"
    trimmed_file = in_fq_base + "_trmd" + in_fq_ext
    porechop_log_path = dirOut_porechop + "pore_" + in_fq_base + ".log"
    
    cmd_porechop = (path_to_porechop + args_porechop + " -i " + in_fq_path +
                    " -o " + dirOut_porechop + trimmed_file + 
                    " -t " + nb_threads)
    
    if TRIM == "porechop":
        TRIM_TIME = t.time()
        Popen_porechop = sub.Popen(cmd_porechop.split(), stdout=sub.PIPE)                
        out_porechop = Popen_porechop.communicate()[0]
        
        with open(porechop_log_path, 'w') as porechop_log:
            porechop_log.write(out_porechop.decode())
            porechop_log.write("TIME TRIM = " + str(t.time()-TRIM_TIME)+"\n\n")
    
    

    # Reads overlapping with Minimap2 followed by chim detection yacrd:
    to_minimap2 = 'minimap2' # If nothing specified, use $PATH Unix var
    if len(list_to_tool) == 1: # NO tool path defined within conf file
        to_minimap2 = list_to_tool[0]
        assert(osp.isfile(to_minimap2))

    args_minimap2_ovlp = "-t" + nb_threads + " -x ava-ont "
    dirOut_yacrd = path_proj + "2_chim/"
    overlapped_file = dirOut_yacrd + in_fq_base + "_ovlp"
    ovlp_paf_path = overlapped_file + ".paf"
    trmd_file_path = dirOut_porechop +  trimmed_file
    flag_ext = "_filtered" + in_fq_ext # Used later and for mapping

    cmd_overlap = (to_minimap2 + "minimap2 " + args_minimap2_ovlp + 
                   trmd_file_path + " " + trmd_file_path)
    
    if CHIM == "yacrd":
        if not osp.isfile(overlapped_file + '.paf'): # 1st need to ovlp
            OVLP_TIME = t.time()
            
            ovlp_log = open(overlapped_file + ".log", 'w')
            ovlp_paf = open(ovlp_paf_path, 'w')
            
            sub.Popen(cmd_overlap.split(), stdout=ovlp_paf, 
                      stderr=ovlp_log).communicate()             
            ovlp_log.write("\n\nTIME OVERLAPPING = " + 
                           str(t.time()-OVLP_TIME) + "\n\n")
            ovlp_log.close()
            ovlp_paf.close()
                
                
        # We can run chimera detection
        YACRD_TIME = t.time()
        path_to_yacrd = path_apps + "miniconda3/envs/metage2019/bin/yacrd"
        log_yacrd_path = dirOut_yacrd + in_fq_base + ".yacrd"
        
        cmd_yacrd = (path_to_yacrd + " -i " + overlapped_file + '.paf' + 
                     " -f " + trmd_file_path) 
        log_yacrd = open(log_yacrd_path, 'w')

        # Start by writting the header:            
        log_yacrd.write("type_of_read" +'\t' + "id_in_mapping_file" +
                        '\t' + "length_of_read" + '\t' +
                        "length_of_gap,begin_pos_of_gap," +
                        "end_pos_of_gap;length_of_gap,dunno,dunno\n")
        Popen_yacrd = sub.Popen(cmd_yacrd.split(), stdout=sub.PIPE)
        log_yacrd.write(Popen_yacrd.communicate()[0].decode())
        # Write elapsed time:
        log_yacrd.write("\n\nTIME YACRD = " + 
                        str(t.time()-YACRD_TIME)+ "\n\n")
        log_yacrd.close()
                     
        # Move the yacrd output fq file to the dirOut_yacrd directory:
        os.rename(dirOut_porechop + osp.splitext(trimmed_file)[0] + 
                  flag_ext, dirOut_yacrd + in_fq_base + flag_ext)           
        
    else:
        pass            
    
   
    # TAXONOMIC DETERMINATION STEP:
    print("Proceeding...")
    
    if DETER == "minimap2":
        dirOut_minimap = outDir
        args_minimap2_map = "-K300M -N25 -t" + nb_threads + " -ax map-ont"
        # args_minimap2_map = "-N25 -t" + nb_threads + " -ax map-ont "
        print("  >> WITH PARAM:", args_minimap2_map)

        root_minimap_outfiles = (dirOut_minimap + in_fq_base + "_to" + 
                                 DB_NAME.capitalize())
        cmd_minimap = " ".join([to_minimap2, args_minimap2_map, 
                                to_ref_db, in_fq_path])
        print() ; print("Core cmd ran:", cmd_minimap)

        START_MINIMAP = t.time()
        log_minimap = open(root_minimap_outfiles + ".log", 'w')
        with open(root_minimap_outfiles + ".sam", 'w') as minimap_SAM:
            sub.Popen(cmd_minimap.split(), stdout=minimap_SAM,
                                           stderr=log_minimap).communicate()
        
        log_minimap.close()
        print("MAPPING FINISHED !")
        
    
    elif DETER == "centri": # Classification using centrifuge:
        dirOut_centri = outDir
        centri_outfile_root = (dirOut_centri + DETER + '_' + in_fq_base + 
                               "_to" + DB_NAME.capitalize())
        
        param_infile = ' -q '
        if in_fq_ext.lstrip('.') in ('fasta', 'fa', 'mfasta'):
            param_infile = ' -f '
        cmd_centri = ("centrifuge -t -p " + nb_threads + param_infile + 
                      in_fq_path + " -x " + to_ref_db +
                      " --report-file " + centri_outfile_root + "_report.tsv " + 
                      "-S " + centri_outfile_root + "_classif.tsv")
        print() ; print("Core cmd ran:", cmd_centri)
        centri_log_path = (dirOut_centri + in_fq_base + "_to" + 
                           DB_NAME.capitalize() + "."  + DETER + "log")
        
        with open(centri_log_path, 'w') as centri_log:
            sub.Popen(cmd_centri.split(),
                      stderr=centri_log).communicate()
        