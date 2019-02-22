#!/usr/bin/env python3


"""
Contains all functions linked with checking of arguments given as input
"""


import os, sys


def check_input_fq(input_fq_path):
    """
    Check if everything's fine with the fq file given as input
    """
    if not os.path.isfile(input_fq_path):
        print("ERROR! Wrong specified input fastq file")
        sys.exit(2)
    else:
        head_input_fq, tail_input_fq = os.path.split(input_fq_path)
        input_fq_base, input_fq_ext = os.path.splitext(tail_input_fq)
        
        if input_fq_ext[1:] not in ("fq", "fastq"):
            err_ext = ("ERROR! Current only fastq file (ext=fq or " +
                       "ext=fastq) are accepted")
            print(err_ext)
            sys.exit(2)
        
        else:
            return (input_fq_path, input_fq_base)


def check_input_nb(input_nb):
    """
    Check if a integer given as argument is valid
    """
    try:
        int(input_nb)
        
    except ValueError:
        print("ERROR! You should give an integer here")
        sys.exit(2)
    
    return input_nb
    

def check_input_taxo_cutoff(input_cutoff, list_taxo):
    """
    Check if a the cutoff for taxonomic level given as argument is valid
    """
    cutoff = input_cutoff.lower()
    if cutoff in list_taxo + ["none"]:
        return input_cutoff
    else:
        print("ERROR! You must enter a valid cutoff for taxonomic level")
        print("Must be in:", list_taxo, '\n')
        sys.exit(2)

