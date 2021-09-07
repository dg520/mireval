#!/usr/bin/env python

from __future__ import division
import os, sys, getopt
import shlex, subprocess
import commands
import numpy
import math


def usage():
    print """
USAGE: ./rfam_chr.py -s <species_name> -o <output_file>

Options:
-s              Species name connected by "_". e.g.: mus_musculus [Compulsory]
-o              Output file name. [Compulsory]
-h              Print this help
--help

"""
    return 1

def par(argv):
    animal=""
    output_loc=""
    try:
        opts, args = getopt.getopt(argv,"s:o:",["sfile=","ofile="])
        #print opts, args
    except getopt.GetoptError:
        usage() # -X : Promyelocyte_dominent_ir_boundary; -Y : Granulocyte_dominent_ir_boundary
        sys.exit(2)
    for opt, arg in opts:
        #print opt,arg
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-s", "--sfile"):
            animal = arg
        elif opt in ("-o", "--ofile"):
            output_loc= arg
        
        else:
            usage()
            sys.exit(2)

    return (animal, output_loc)


cwd = os.getcwd()
root="/".join(cwd.split("/")[:-1])

species,output_dir =par(sys.argv[1:])

if species !="" and output_dir!="":
    if species !="rfam" and species !="rfam":
        spp=species.split("_")[0]+" "+species.split("_")[1]
        spp1="chromosome"
        run="for i in "+root+"/mireval/ref/rfam/*.gff3;do echo $i;head -n 10 $i | grep -m 1 -i \""+spp+"\";done >"+"toto"
        #print run
        sta,out=commands.getstatusoutput(run)  
        run="less toto | grep -B 1 -i \""+spp1+"\""
        sta,out=commands.getstatusoutput(run)
        p=out.split("\n")
        rf=[]
        if "chromosome " in out:
            i=0
            for x in p:
                if x!="--":
                    if i%2==0:
                        tmp=x
                    if i%2==1:
                        #print x
                        rf.append("chr"+x.split("chromosome ")[1].split(",")[0]+"\t"+tmp+"\n")
                    i+=1
        if "chromosome, " in out:
            i=0
            for x in p:
                if x!="--":
                    if i%2==0:
                        tmp=x
                    if i%2==1:
                        #print x
                        rf.append(x.split("chromosome, ")[1].split(",")[0]+"\t"+tmp+"\n")
                    i+=1
        f=open(output_dir,"w")
        f.writelines(rf)
        f.close()

        run="rm toto"
        sta,out=commands.getstatusoutput(run)

else:
    usage()
    if species=="":
        print "Error: Need a species name to extract."
    if output_dir=="":
        print "Error: Need an output file name."
