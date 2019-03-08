#!/usr/bin/env python3

"""
Take a mfasta, containing the sequences associated with each copy of the 16S,
of a given organism
"""

import os, sys, re
import os.path as osp
import itertools as ittls
import subprocess as sub
from Bio import SeqIO
import pandas as pd


def align_two_seq(to_fasta_copy, to_main_fasta):
    cmd_needle = ("needle -asequence " + to_fasta_copy + " -bsequence " + 
                  to_main_fasta + " -out stdout -gapopen 10 -gapextend 0.5")
    salut = sub.Popen(cmd_needle.split(), stdout=sub.PIPE,
                      stderr=sub.DEVNULL).communicate()[0].decode()

    my_regex = re.compile("(?:# Identity:\s+)(.*)")
    matches = my_regex.search(salut)
    id_perc = float(matches.group(1).split()[-1].rstrip('%)').lstrip('('))

    return id_perc


def print_upper_arr(df_arr):
    """
    Display only upper matrix
    """
    print("", end='\t') # Avoid auto newline printing
    for col_name in df_arr.columns:
        print(col_name, end='\t')
    print()
    df_arr = df_arr.reset_index()
    nb_col = len(df_arr.columns)
    for i in range(nb_col-1):
        print(df_arr.loc[i, 'index'], end='\t')
        for j in range(nb_col-1):
            if  j >= i:
                print(df_arr.loc[i, "copy_"+str(j+1)], end='\t')
            else:
                print("", end='\t')
        print()         




# MAIN:
if __name__ == '__main__':
    # to_16S_fa = "/home/sheldon/SEGO/V1_SC/16S/Bacillus_subtilis_16S.fasta"
    to_16S_fa = sys.argv[1]
    filename = osp.basename(to_16S_fa)
    base_in_fasta = filename.split('.')[0]
    inter_mode = True # Mode

    name_tmp_folder = "tmp_needle"
    if osp.isdir(name_tmp_folder):
        os.system("rm -r " + name_tmp_folder)
    os.mkdir(name_tmp_folder)

    list_mfasta_rec = list(SeqIO.parse(to_16S_fa, 'fasta'))
    nb_copies_16S = len(list_mfasta_rec)
    list_names = []
    iterator = ittls.combinations_with_replacement(range(1, nb_copies_16S+1), 2)
    df_id_percents = pd.DataFrame(data=None)

    # Write all tmp fasta files (1 for each copy of 16S):
    for rec_fasta in list_mfasta_rec:
        list_names.append(rec_fasta.name)
        out_fasta = osp.join(name_tmp_folder, rec_fasta.name)
        SeqIO.write(rec_fasta, out_fasta, "fasta-2line")

    for row, col in iterator:
        name_row = list_names[row-1]
        out_fasta_row = osp.join(name_tmp_folder, name_row)

        name_col = list_names[col-1]
        out_fasta_col = osp.join(name_tmp_folder, name_col)

        if inter_mode:
            idx_row = name_row[0:5] + "_" + name_row.split('_')[-1]
            idx_col = name_col[0:5] + "_" + name_col.split('_')[-1]
        else: # Intra-sp mode
            idx_row = "copy_" + name_row.split('_')[-1]
            idx_col = "copy_" + name_col.split('_')[-1]

        print("Alignment of:", name_row, "and", name_col, "...")
        id_percent = align_two_seq(out_fasta_row, out_fasta_col)

        if idx_row == idx_col: # Check autosimil == 100%
            assert(int(id_percent) == 100)

        df_id_percents.loc[idx_row, idx_col] = id_percent

    os.system("rm -r " + name_tmp_folder) # Cleaning

    # Printing of results:
    print("\n\n")
    if inter_mode:
        print("", end='\t')
        [print(row_name, end=' ') for row_name in df_id_percents.index]
        print()
        import numpy as np
        np.set_printoptions(linewidth=100) # To avoid auto-wrapping at 75 char
        for row_name in df_id_percents.index:
            print(row_name, end='  ') 
            print(np.asarray(df_id_percents.loc[row_name, :]))
    else:
        print(base_in_fasta, ":")
        print_upper_arr(df_id_percents)
    print()