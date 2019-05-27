#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: add_gene_annotation2deseq2_results.py
- CONTACT: Gaurav Jain(gaurav.jain@dzne.edu)
***********************************************
"""
print (__doc__)

# Built in modules
import argparse
import os.path
import sys

# 3rd party modules
import textwrap
import re
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mp
#mp.use('Agg') # to use matplotlib without X11
import matplotlib.pyplot as plt
import subprocess
import binascii as bi
import scipy.stats as stats
from collections import *
from numpy import nanmean

# for looping files in a dir
import glob

# user defined modules
from gjainLIB import *      # import all the functions from the Gaurav`s python library

### for color scale
from  matplotlib import colors
from itertools import cycle, islice # barplot colors

################ USER CONFIGURATION ###################
np.set_printoptions(precision=6)
#######################################################

def main():
    # Get input options
    args = check_options()

    # Store the variables
    annotation_file  = args.annotation_file
    input_file       = args.input_file
    output_file      = args.output_file

    
    # Get output file name
    if not output_file:
        output_file = "{0}_geneSymbols.{1}".format(get_file_info(input_file)[3], get_file_info(input_file)[2])
    ofile = open(output_file, 'w')

    # Get gene annotations
    gannotation = get_gene_annotation_dict(annotation_file)

    print "- Reading and add annotation to file..."
    with open(input_file, 'rU') as f:
        # Get the header
        header  = f.readline().strip()
        hrow    = re.split("\t", header)
        hrow.insert(0, "GeneName")
        hline = "\t".join(hrow)
        ofile.write("{0}\n".format(hline))

        for line in (l.strip() for l in f if not l.startswith('#')):
            # head output/magda/FC_day_night/mitoutliers/01.2_DEseq/results_DEseq2/CSnCA1_over_CnCA1_DE_RESULTS.txt | cut -f1-3
            # feature	baseMean	log2FoldChange
            # ENSMUSG00000084164	1.69825668905681	0.944893830777482
            # ENSMUSG00000063830	2.15921481366273	0.911525811293581
            # .....

            # Split the line
            row    = line.split("\t") 
            if not row[0] in gannotation:
                name = ["no_annotation", row[0]]
            else:
                name = re.split("\|", gannotation[row[0]])         # Kitl|ENSMUSG00000019966|ENSMUST00000105283
            
            row.insert(0,name[0]) 
            line  = "\t".join(row)
            ofile.write("{0}\n".format(line))
    ofile.close()
    
    # Convert to xls file
    print "\t- Converting to excel..."
    os.system("ssconvert {0} {1}.xls".format(output_file, get_file_info(output_file)[3]))

    print "\n- Your text  output file is: {0}".format(output_file)
    print "\n- Your excel output file is: {0}.xls".format(get_file_info(output_file)[3])
            
################ USER DEFINED FUNCTIONS ###################
def get_gene_annotation_dict(input_file):
    ''' Get the annotation in a dictionary '''

    # Get the information in a dictionary
    d = defaultdict(int)
    with open(input_file,'rU') as fg:
        # head input/annotation/Mus_musculus_Ensembl_GRCm38.genes.bed 
        # chr10	100015630	100016020	Kitl|ENSMUSG00000019966|ENSMUST00000105283	0	+
        # chr10	100015630	100016035	Kitl|ENSMUSG00000019966|ENSMUST00000105283	0	+
        # chr10	100015630	100100413	Kitl|ENSMUSG00000019966|ENSMUST00000105283	0	+
        # chr10	100015826	100016020	Kitl|ENSMUSG00000019966|ENSMUST00000130190	0	+
        # ....
        
        # Loop through rest of the file
        for line in useful_lines(fg):
            name            = re.split('\t', line)[3]
            geneid          = re.split('\|', name)[1]
            transcriptid    = re.split('\|', name)[2]
            d[geneid]       = name

    return d

def print_help():
    ''' Print system help '''
    print >> sys.stderr, "\n ----------------- HELP ------------------\n", parser.print_help(), "\n"

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/add_gene_annotation2deseq2_results.py -af=input/annotation/Mus_musculus_Ensembl_GRCm38.genes.bed -if=input/exomePeaks/1week_sig_diff_peak.txt
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-if", metavar='--infile', help="input  file", dest="input_file", type=str, required=True)
    parser.add_argument("-af", metavar='--anfile', help="input  annotation file", dest="annotation_file", type=str, required=True)
    parser.add_argument("-of", metavar='--otfile', help="output file", dest="output_file", type=str)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().input_file:
        logdir = "{0}/logs".format(get_file_info(parser.parse_args().input_file)[0])
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir, get_file_info(parser.parse_args().input_file)[1])
    else:
        logdir  = "{0}/logs".format(os.getcwd())
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir,get_file_info(sys.argv[0])[1])

    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args
    
if __name__=="__main__":
      main()


