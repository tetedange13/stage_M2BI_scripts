#!/usr/bin/env python3


"""
Contains all functions linked with checking of arguments given as input
"""


import os, sys


def infile(infile_path, list_acceptable_ext):
    """
    Check if everything's fine with the fq or fa file given as input
    """
    if not os.path.isfile(infile_path):
        print("ERROR! Specified input file does not exit")
        print() ; sys.exit(2)
    else:
        head_infile, tail_infile = os.path.split(infile_path)
        infile_base, infile_ext = os.path.splitext(tail_infile)
        
        if infile_ext[1:] not in list_acceptable_ext:
            err_ext = ("ERROR! The given file as the wrong extension. " +
                       "Should be in:", list_acceptable_ext)
            print(err_ext)
            sys.exit(2)
        
        else:
            return (infile_path, infile_base, head_infile, tail_infile, 
                    infile_ext)
            
            
def input_nb(input_nb, param_name=""):
    """
    Check if a integer given as argument is valid
    """
    try:
        int(input_nb)
        
    except ValueError:
        print("ERROR! You should give an integer for the parameter:", 
              param_name)
        print() ; sys.exit(2)
    
    return input_nb


def bool_type(rep):
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
        print() ; sys.exit(2)


def acceptable_str(input_str, list_acceptable):
    """
    Check if a string given by the user is valid, i.e. belongs to the list of
    acceptable values
    """
    in_str = input_str.lower()
    if in_str in list_acceptable:
        return in_str
    else:
        print("ERROR! String given as parameter not valid. Should be in:",
              list_acceptable)
        print() ; sys.exit(2)


def list_config(list_config, name_param):
    """
    Take a list of config arguments and check its size (supposed to be of 
    size 1), to decided whether to return the path (within the list) or send
    an error
    """
    if len(list_config) != 1:
        if not list_config: # Empty list ==> NO config found for given params:
            print("ERROR: NO config found with given args for '{}' param".format(name_param))
        else: # len() > 1 ==> More than 1 possible config
            print("ERROR: Several config possible for '{}' param".format(name_param))
        print() ; sys.exit(2)
    
    return list_config[0]  
