#!/usr/bin/env python3

"""
Long-reads (LR) classification pipeline

Usage:
  1-assign_pipeline.py (-i <inputFqFile>) (-c <confFile>) (-P <prog>) (-d <db>) [-t <threads>] [-o <outDir>]
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inputFqFile=input_fq  input fastq file
  -c --confFile=conf_file    path to the configuration file
  -d --db=database           name of the database to use
  -P --prog=program          tool to use for taxonomic determination step
  -t --threads=nb_threads    number of threads for parallelisation [default: 10]
  -o --outDir=out_dir        output directory for taxo file (SAM or CSV) [default: ./]
  
Arguments:
    configuration file: e.g. '/projets/metage_ONT_2019/stage_M2BI_scripts/pipeline.conf'
    determination: 'minimap2', 'centri'
    database: 'SILVA', 'rrn', 'pCompressed', 'zymo', 'ncbi16S' or 'newLot_Zymo'
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

    # Check if db_name does not contain any underscore (problematic for 2-prim_analysis.py):
    DB_NAME = ARGS["--db"]
    if DB_NAME.find('_') != -1:
        print("ERROR: Given database name contains (at least) an underscore")
        camel_db = ''.join(x for x in DB_NAME.title() 
                                   if x != '_')
        print("(you should replace it with camelcase: {})".format(camel_db))
        print() ; sys.exit()

    PROG = check.acceptable_str(ARGS["--prog"], ["minimap2", "centri"])
    nb_threads = check.input_nb(ARGS["--threads"], "'-t number of threads'")
    outDir = ARGS["--outDir"]
    if not osp.isdir(outDir):
        print("ERROR OutDir: {} does NOT exit !\n".format(outDir)) ; sys.exit()
    to_conf_file = ARGS['--confFile']
    if not osp.isfile(to_conf_file):
        print("ERROR ConFile: {} does NOT exit !\n".format(conf_file)) ; sys.exit()

    print()
    print("Tool: {} | Against: {}".format(PROG, DB_NAME))
    print("Output directory:", outDir)


    # Extracting parameters from 'pipeline.conf' file:
    conf_csv = pd.read_csv(to_conf_file, sep=';', comment='#')
    # Lowercase conversion, to make it more flexible:
    conf_csv.tool.str.lower() ; conf_csv.type_param.str.lower()
    params_csv = conf_csv[conf_csv.tool == PROG].drop(['tool'], axis='columns')

    cond_to_tool = params_csv.type_param == 'to_exe'
    list_to_tool = params_csv[cond_to_tool].supplField_1.values
    if len(list_to_tool) > 1:
        print("ERROR: Several possible params for 'to_exe'") ; sys.exit()

    cond_to_ref_db = ((params_csv.type_param == 'to_ref') & 
                      (params_csv.supplField_1.str.lower() == DB_NAME.lower()))
    list_to_ref_db = params_csv[cond_to_ref_db].supplField_2.values
    to_ref_db = check.list_config(list_to_ref_db, 'to_ref')
    if PROG != 'centri':
        examined_file = to_ref_db
    else:
        examined_file = to_ref_db + '.1.cf'

    if not  osp.isfile(examined_file):
        print("ERROR: Ref file (for DB) provided does NOT exist:", 
              examined_file)
        print() ; sys.exit(2)
    
    print()
    print("Deducted from conf file:")
    print("To reference database:", to_ref_db)
    

    # Determination of the path towards Minimap2 executable:
    to_minimap2 = 'minimap2' # If nothing specified, use $PATH Unix var
    if len(list_to_tool) == 1: # Else use at the given path
        to_minimap2 = list_to_tool[0]
        assert(osp.isfile(to_minimap2))

   
    # TAXONOMIC DETERMINATION STEP:
    print("Proceeding...")
    
    if PROG == "minimap2":
        dirOut_minimap = outDir
        # The 'MD' tag is requiered for some tools
        args_minimap2_map = "-K300M -N25 --MD -t" + nb_threads + " -ax map-ont"
        # args_minimap2_map = "-N25 -t" + nb_threads + " -ax map-ont "
        print("  >> WITH PARAM:", args_minimap2_map)

        root_minimap_outfiles = (dirOut_minimap + in_fq_base + "_to" + 
                                 DB_NAME[0].upper() + DB_NAME[1:])
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
        
    
    elif PROG == "centri": # Classification using centrifuge:
        dirOut_centri = outDir
        centri_outfile_root = (dirOut_centri + PROG + '_' + in_fq_base + 
                               "_to" + DB_NAME[0].upper() + DB_NAME[1:])
        
        param_infile = ' -q '
        if in_fq_ext.lstrip('.') in ('fasta', 'fa', 'mfasta'):
            param_infile = ' -f '
        cmd_centri = ("centrifuge -t -p " + nb_threads + param_infile + 
                      in_fq_path + " -x " + to_ref_db +
                      " --report-file " + centri_outfile_root + "_report.tsv " + 
                      "-S " + centri_outfile_root + "_classif.tsv")
        print() ; print("Core cmd ran:", cmd_centri)
        centri_log_path = centri_outfile_root + ".log"
        
        with open(centri_log_path, 'w') as centri_log:
            sub.Popen(cmd_centri.split(),
                      stderr=centri_log).communicate()
        
