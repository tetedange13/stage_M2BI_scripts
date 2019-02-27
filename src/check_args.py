#!/usr/bin/env python3


"""
Contains all functions linked with checking of arguments given as input
"""


import os, sys


def check_infile(infile_path, list_acceptable_ext):
    """
    Check if everything's fine with the fq or fa file given as input
    """
    if not os.path.isfile(infile_path):
        print("ERROR! Wrong specified input fastq file")
        sys.exit(2)
    else:
        head_infile, tail_infile = os.path.split(infile_path)
        infile_base, infile_ext = os.path.splitext(tail_infile)
        
        if infile_ext[1:] not in list_acceptable_ext:
            err_ext = ("ERROR! The given file as the wrong extension. " +
                       "Should be in:" + list_acceptable_ext)
            print(err_ext)
            sys.exit(2)
        
        else:
            return (infile_path, infile_base, head_infile, tail_infile, 
                    infile_ext)
            
            
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


def check_bool_type(rep):
    """
    Check if the answer is boolean type.
    Args:
        rep: answer given by the user (str)
    Returns:
        The answer casted into a boolean
    """
    rep = rep.lower()
    if rep in ("true", "t", "yes", "y"):
        return True
    elif rep in ("false", "f", "no", "n"):
        return False
    else:
        print("Enter a boolean type (True, T, False, F) !")
        sys.exit(2)
        
        
def check_soft(trim, list_softs):
    """
    Check if a suitable name of soft has been given as argument
    """
    trim = trim.lower()
    if trim in ["no"] + list_softs:
        return trim
    else:
        print("Enter a valid soft name or 'no' (default) !")
        sys.exit(2)
        
