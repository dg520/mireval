from __future__ import division
import os, sys, getopt
import shlex, subprocess
import commands
import numpy
import math
from Bio import SeqIO


def usage():
    print """
USAGE: ./countChromSizes.py -i <Genome_FATSTA_file> -o <output_file_name>

Options:
-i              Single FASTA file for the whole Genome. [Compulsory]
-o              Output file name. [Compulsory]
-h              Print this help
--help

"""
    return 1


def par(argv):
    input_loc=""
    output_loc=""
    try:
        opts, args = getopt.getopt(argv,"i:o:",["ifile=","ofile="])
        #print opts, args
    except getopt.GetoptError:
        usage() # -X : Promyelocyte_dominent_ir_boundary; -Y : Granulocyte_dominent_ir_boundary
        sys.exit(2)
    for opt, arg in opts:
        #print opt,arg
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_loc = arg
        elif opt in ("-o", "--ofile"):
            output_loc= arg
        
        else:
            usage()
            sys.exit(2)

    return (input_loc, output_loc)


input_dir,output_dir =par(sys.argv[1:])
if input_dir !="" and output_dir!="":
    handle = open(input_dir, "rU")
    f=open(output_dir,"w")
    for record in SeqIO.parse(handle, "fasta"):
        f.write(record.id+"\t"+str(len(record.seq))+"\n")
    handle.close()
    f.close()


else:
    usage()
    if input_dir =="":
        print "Error: Need an input FASTA file for the genome."
    if output_dir =="":
        print "Error: You should indicate an output file name."
