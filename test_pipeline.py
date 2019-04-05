#!/usr/bin/env python3

"""
Long-reads (LR) classification pipeline

Usage:
  main.py (-i <inputFqFile>) [-d <db>] [-T <trim>] [-C <chim>] [-D <deter>] [-t <threads>]
  
Options:
  -h --help                  help
  --version                  version of the script
  -i --inputFqFile=input_fq  input fastq file
  -d --db=database           name of the database to use [default: Zymo]
  -T --trim=trimming         use Porechop to trim reads [default: no]
  -C --chim=chimera          chimera detection using yacrd [default: no]
  -D --deter=determination   taxonomic determination step [default: no]
  -t --threads=nb_threads     number of threads for parallelisation [default: 10]
  
Arguments:
    trimming: 'porechop' or 'no' (default)
    chimera: 'yacrd' or 'no' (default)
    determination: 'minimap2', 'margin', 'centri', or 'no' (default)
    database: 'GG', 'SILVA', rrn' or 'Zymo'
"""

import os
import sys
import subprocess as sub
import time as t
import os.path as osp
from docopt import docopt
import src.check_args as check
        
        

# MAIN:
if __name__ == "__main__":
    START_TIME = t.time()
    
    # Get arguments:
    ARGS = docopt(__doc__, version='0.1')
    tuple_check_fq = check.infile(ARGS["--inputFqFile"], ('fq', 'fastq'))
    in_fq_path, in_fq_base, _, tail_input_fq, in_fq_ext = tuple_check_fq
    DB_NAME = check.acceptable_str(ARGS["--db"], 
                                   ["rrn", "zymo", "silva", "gg"])
    TRIM = check.acceptable_str(ARGS["--trim"], ["porechop", "no"])
    CHIM = check.acceptable_str(ARGS["--chim"], ["yacrd", "no"])
    DETER = check.acceptable_str(ARGS["--deter"], 
                                 ["minimap2", "margin", "centri", "no"])
    nb_threads = check.input_nb(ARGS["--threads"], "'-t number of threads'")   

                
    #Common variables/params:
    # To databases directory:
    to_dbs = "/mnt/72fc12ed-f59b-4e3a-8bc4-8dcd474ba56f/metage_ONT_2019/"
    
    path_apps = "/home/sheldon/Applications/"
    path_proj = "/projets/metage_ONT_2019/"
    #print("Processing " + tail_input_fq, '\n')
    
    
    # Adaptators trimming using Porechop:
    path_to_porechop = "/home/sheldon/SEGO/Porechop/porechop-runner.py"
    #args_porechop = " --min_split_read_size 100 --extra_middle_trim_bad_side 10"
    #args_porechop = " --middle_threshold 90"
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
        
    else:
        pass
        # print('\n', cmd_porechop)
        # print(porechop_log_path, '\n')

    #print("Trimming finished !")
    
    

    # Reads overlapping with Minimap2 followed by chim detection yacrd:
    to_marginAlign = path_apps + "marginAlign-23jan19/" # Minimap2 inside
    to_minimap2 = to_marginAlign + "submodules/minimap2/"
    #path_to_paftools = to_minimap2 + "misc/"
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
        # print(cmd_yacrd);sys.exit()
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
        # print('\n', cmd_overlap)
        # print(overlapped_file, '\n')               
    
   
    # TAXONOMIC DETERMINATION STEP:
    print("TAXONOMIC DETERMINATION WITH:", DETER + "...")
    file_to_map = dirOut_yacrd + in_fq_base + flag_ext
    
    if DETER == "margin":
        # Mapping using MarginAlign+minimap2:
        dirOut_margin = path_proj + "3-deter_margin/"
        args_margin = " --minimap2 --maxThreads=" + nb_threads + " "
        base_margin = dirOut_margin + in_fq_base
        
        cmd_margin = (to_marginAlign + "marginAlign" + args_margin +
                           file_to_map + " " + ref_fa_path + " " + 
                           base_margin + "_toZymo.sam" +
                           " --em --iterations=30 " + 
                           "--maxAlignmentLengthToSample=50000000 " +
                           "--outputModel " + 
                           base_margin + ".hmm" + " --jobTree=" +
                           base_margin + "_jobTree " + "--trials=5")
        margin_log_path = dirOut_margin + in_fq_base + "_margin.log"
        
        #print('\n', cmd_margin, '\n') ; sys.exit()
        #print("Mapping with MarginAlign+minimap2 in progress...")
        MARGIN_TIME = t.time()
        with open(margin_log_path, 'w') as margin_log:
            sub.Popen(cmd_margin.split(), stderr=margin_log).communicate()                
        #out_margin, err_margin = Popen_margin.communicate()
        
        #margin_log = open(margin_log_path, 'w')
        #if out_margin:
        #    margin_log.write(out_margin)
        #if err_margin:
        #    margin_log.write(err_margin)
        
            margin_log.write("\n\nTIME MARGIN = " + str(t.time()-MARGIN_TIME) + 
                             "\n\n")       
        #margin_log.close()        
    
    #print("Mapping with MarginAlign+minimap2 finished !")
    
    elif DETER == "minimap2":
        # dirOut_minimap = path_proj + "3-deter_minimap2/"
        dirOut_minimap = path_proj + "3-deter_minimap2_second/"

        args_minimap2_map = "-N5000 -t" + nb_threads + " -x map-ont "
        # if DB_NAME == "gg":
        #     args_minimap2_map = ("-K100M --secondary=no -t" + nb_threads + 
        #                          " -x map-ont ")

        to_minimap_idxes = to_dbs + "minimap2_idxes/"
        root_minimap_outfiles = (dirOut_minimap + in_fq_base + "_to" + 
                                 DB_NAME.capitalize())
        cmd_minimap = (to_minimap2 + "minimap2 " + args_minimap2_map + "-a " +
                       to_minimap_idxes + DB_NAME + ".mmi " + 
                       in_fq_path)
        # print(cmd_minimap) ; sys.exit()
        START_MINIMAP = t.time()
        log_minimap = open(root_minimap_outfiles + ".log", 'w')
        with open(root_minimap_outfiles + ".sam", 'w') as minimap_SAM:
            sub.Popen(cmd_minimap.split(), stdout=minimap_SAM,
                                           stderr=log_minimap).communicate()
        
        #log_minimap.write("\n\nRUNTIME MAPPING WITH MINIMAP2: " + 
        #                  str(t.time() - START_MINIMAP)) # --> MARCHE PAS
        log_minimap.close()
        print("MAPPING FINISHED !")
        sys.exit()
        
    
    elif DETER == "centri": # Classification using centrifuge:
        # CMD COMPILATION CENTRIFUGE:
#/usr/bin/g++  -O3 -m64 -msse2 -funroll-loops -g3 -std=c++11 -DCOMPILER_OPTIONS="\"-O3 -m64 -msse2 -funroll-loops -g3 -std=c++11 -DPOPCNT_CAPABILITY\"" -DPOPCNT_CAPABILITY -fno-strict-aliasing -DCENTRIFUGE_VERSION="\"1.0.4\"" -DBUILD_HOST="\"`hostname`\"" -DBUILD_TIME="\"`date`\"" -DCOMPILER_VERSION="\"`/usr/bin/g++  -v 2>&1 | tail -1`\"" -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE   -DBOWTIE_MM   -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX -DNDEBUG -Wall  -I third_party  -o centrifuge-class centrifuge.cpp ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp edit.cpp bt2_idx.cpp reference.cpp ds.cpp limit.cpp random_source.cpp tinythread.cpp qual.cpp pat.cpp read_qseq.cpp ref_coord.cpp mask.cpp pe.cpp aligner_seed_policy.cpp scoring.cpp presets.cpp simple_func.cpp random_util.cpp outq.cpp centrifuge_main.cpp -lpthread
        
        to_Centri_idxes = to_dbs + "Centri_idxes/"
        dirOut_centri = path_proj + "3-deter_centri/"
        centri_classif_path = (dirOut_centri + in_fq_base + "_to" +  
                               DB_NAME.capitalize() + "_classif.tsv")
        
        cmd_centri = ("centrifuge -t -p " + nb_threads + " -q " + 
                      in_fq_path + " -x " + to_Centri_idxes + DB_NAME +
                      " --report-file " + dirOut_centri + in_fq_base + "_to" + 
                      DB_NAME.capitalize() + "_report.tsv -S " +
                      centri_classif_path)
        # print(cmd_centri);sys.exit()
        
        centri_log_path = (dirOut_centri + in_fq_base + "_to" + 
                           DB_NAME.capitalize() + "."  + DETER + "log")
        
        # with open(centri_classif_path, 'w') as classif_centri, \
        #      open(centri_log_path, 'w') as centri_log:
        with open(centri_log_path, 'w') as centri_log:
            sub.Popen(cmd_centri.split(),
                      stderr=centri_log).communicate()
    
    
    sys.exit()

    # De novo clustering using CARNAC-LR:
    path_to_carnac = path_apps + "CARNAC-LR_git93dd640/"
    dirOut_carnac = "./CARNAC/"
    carnac_filout = dirOut_carnac + in_fq_base # Prefix for CARNAC output
    cmd_carnac_convert = ("python3 " + path_to_carnac + 
                          "scripts/paf_to_CARNAC.py " + ovlp_paf_path + " " +
                          trmd_file_path + " input_CARNAC.txt")                   
    cmd_carnac = (path_to_carnac + "CARNAC-LR -f input_CARNAC.txt -o " + 
                  carnac_filout + "_CARNAC.txt -t " + 
                  nb_threads)
    
    proceed = False
    if proceed:
        CARNAC_TIME = t.time()
        os.system(cmd_carnac_convert)
        
        with open(carnac_filout + "_CARNAC.log", 'w') as carnac_log:
            sub.Popen(cmd_carnac.split(), stdout=carnac_log).communicate()
            carnac_log.write("\n\nTIME CARNAC = " + str(t.time()-CARNAC_TIME) + 
                             "\n\n")
         
        os.remove("input_CARNAC.txt") # input_CARNAC is temporary
        # Moving CARNAC "*_metrics.txt" files:
        os.rename("clusters_metrics.txt", 
                  dirOut_carnac + in_fq_base + "_clust_metrics.txt")
        os.rename("nodes_metrics.txt", 
                  dirOut_carnac + in_fq_base + "_nodes_metrics.txt")
    
    else:
        print(cmd_carnac_convert)
        print(cmd_carnac)
    
    sys.exit()
        
    # Mapping using NGMLR:
    NGMLR_TIME =  t.time()
    path_to_ngmlr = (path_apps + "ngmlr-0.2.7/bin/ngmlr-0.2.7/"+
                     "ngmlr")
    cmd_ngmlr = (path_to_ngmlr + " -t 4 -r reference.fasta -q reads.fastq " +
                 "-o test.sam -x ont")
    print("Mapping with NGMLR in progress...")
    print('\n', cmd_ngmlr, '\n')
    print("Mapping with NGMLR finished !")
    print("TIME NGMLR =", t.time()-NGMLR_TIME, '\n')
    
    
    
    # Mapping with LordFAST:*
    path_to_lordfast = path_apps + "lordfast-0.0.10/lordfast"
    #cmd_lordfast = (path_to_lordfast + "--search " + refgen.fa " --seq " +
    #                reads.fastq + " > map.sam")
        