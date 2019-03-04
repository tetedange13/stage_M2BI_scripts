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
import numpy as np
import pandas as pd

def needle_align(seq_a, seq_b):
#def align_two():
    cmd_needle = ("needle -asequence " + seq_a + " -bsequence " + seq_b)
    print(cmd_needle)
    return 0.8

def align_two_seq(to_fasta_copy, to_main_fasta):
    cmd_needle = ("needle -asequence " + to_fasta_copy + " -bsequence " + 
                  to_main_fasta + " -out stdout -gapopen 10 -gapextend 0.5")
    salut = sub.Popen(cmd_needle.split(), 
                      stdout=sub.PIPE, 
                      stderr=sub.DEVNULL).communicate()[0].decode()
    # my_regex = re.compile("(?:# 1: )(.*)(?:\n# 2: )(.*)(?:(\n.*){6})(.*)")
    my_regex = re.compile("(?:# Identity:\s+)(.*)")
    # all_matches = my_regex.findall(salut)
    matches = my_regex.search(salut)
    # print(toto.group(1));sys.exit()
    # seq_row, seq_col = matches.group(1), matches.group(2)
    # id_percent = float(matches.group(3).split()[-1].rstrip('%)').lstrip('('))
    id_perc = float(matches.group(1).split()[-1].rstrip('%)').lstrip('('))
    # #to_return = [match[0], match[for match in ]
    # return (seq_row, seq_col, id_percent)
    # return all_matches
    return id_perc


# MAIN:
if __name__ == '__main__':
    to_16S_fa = "/home/sheldon/SEGO/V1_SC/16S/Bacillus_subtilis_16S.fasta"
    filename = osp.basename(to_16S_fa)
    base_in_fasta = filename.split('.')[0]
    # print(base_in_fasta);sys.exit()

    name_tmp_folder = "tmp_needle"
    if osp.isdir(name_tmp_folder):
        os.system("rm -r " + name_tmp_folder)
    os.mkdir(name_tmp_folder)

    list_mfasta_rec = list(SeqIO.parse(to_16S_fa, 'fasta'))
    nb_copies_16S = len(list_mfasta_rec)
    list_names = []
    toto = ittls.combinations_with_replacement(range(1, nb_copies_16S+1), 2)
    df_id_percents = pd.DataFrame(data=None)

    # Write all tmp fasta files (1 for each copy of 16S):
    for rec_fasta in list_mfasta_rec:
        list_names.append(rec_fasta.name)
        out_fasta = osp.join(name_tmp_folder, rec_fasta.name)
        SeqIO.write(rec_fasta, out_fasta, "fasta-2line")

    # sys.exit()
    for row, col in toto:
        name_row = list_names[row-1]
        out_fasta_row = osp.join(name_tmp_folder, name_row)

        name_col = list_names[col-1]
        out_fasta_col = osp.join(name_tmp_folder, name_col)

        idx_row = "copy_" + name_row.split('_')[-1]
        idx_col = "copy_" + name_col.split('_')[-1]

        # if idx_row == idx_col: # No need to align
        if False:
            id_percent = 100.0
        else:
            print("Alignment of:", name_row, "and", name_col, "...")
            id_percent = align_two_seq(out_fasta_row, out_fasta_col)

        df_id_percents.loc[idx_row, idx_col] = id_percent
        # break
        
        # os.remove(out_fasta_row)
        # if osp.isfile(out_fasta_col):
        #     os.remove(out_fasta_file)

        # break
    print("\n\n", df_id_percents)