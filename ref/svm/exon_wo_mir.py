from __future__ import division
import shlex, subprocess
import commands
import cgi, cgitb
from operator import itemgetter
from Bio import Motif
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.Blast import NCBIXML
import re
import random
from random import randint


def neighbour_mir(root_name,chromosome,chromosome0,start_cor,end_cor,strand,species):
    files=open(root_name+"/web/genomes/"+species+"/mir/mir.gff3","r")
    mir_list1=files.readlines()
    files.close()
    mir_list="\n".join(mir_list1).split("#\n\n")[-1].split("\n\n")

    caught=0
    for premir in mir_list:
        if ("miRNA_primary_transcript" in premir)==False:
            mirchr=premir.split("\t")[0]
            if mirchr[0:3]!="chr":
                mirchr="chr"+mirchr
            if mirchr[-2:]==".1":
                mirchr=mirchr[:-2]
            new_chromosome="NA"
            if "_" in chromosome:
                new_chromosome="chr"+chromosome.split("_")[1]
            mirstart=int(premir.split("\t")[3])
            mirend=int(premir.split("\t")[4])
            mirname=premir.split("\t")[-1].split("Name=")[-1].split(";")[0][:-1]
            miracc=premir.split("\t")[-1].split("accession_number=")[-1].split(";")[0]
            mirst=premir.split("\t")[6]
            if (mirchr==chromosome or mirchr==new_chromosome) and mirst==strand:
                if max(start_cor,mirstart)<min(end_cor,mirend):
                    caught=1    
    return caught



def overlap_finder(chromosome,start_cor,end_cor,strand,list1):
    caught=0
    c_t=list1.split("\t")
    c_chr=c_t[0][1:]
    c_start=int(c_t[1])
    c_end=int(c_t[2])
    c_st=c_t[3]
    if c_chr==chromosome and c_st==strand:
        if max(start_cor,c_start)<min(end_cor,c_end):
            caught=1
    return caught
    


root="/home/dadi/projects"
f=open("exonseq83fold.txt")       # a folding file
d=f.readlines()
f.close()

n=0
d1=[]
i_list=[]
while n<len(d):
    i1=randint(0,int(len(d)/3))
    if (i1 in i_list)==False:
        i=i1*3
        t=d[i].split("\t")
        chr=t[0][1:]
        start=int(t[1])
        end=int(t[2])
        st=t[3]
        if neighbour_mir(root,chr,chr,start,end,st,"mus_musculus")>0:
            print d[i]        
        else:
            seq=d[i+1][:-1]
            fold=d[i+2].split(" ")[0]
            fol=d[i+2].split(" ")[0]
            mfe=d[i+2].split(" ")[-1].split(")")[0].split("(")[-1]
            num=0
            for x in fold:
                if x=="(":
                    num+=1
            re1=re.compile("\(\.*\)")
            loop=re1.findall(fol)

            if num>=18 and float(mfe)<-15:
                if d1==[]:
                    d1.append(d[i])
                    d1.append(d[i+1])
                    d1.append(d[i+2])
                else:
                    hit=0
                    a=0
                    while a<len(d1):
                        hit=hit+ overlap_finder(chr,start,end,st,d1[a])
                        a+=3
                    if hit==0:
                        d1.append(d[i])
                        d1.append(d[i+1])
                        d1.append(d[i+2])
    #i+=3
    #print len(d1)
    i_list.append(i1)
    if len(d1)==21279:
        break


f=open("exon_mirlike_mm10_nonredundency_nonmir.txt","w")
f.writelines(d1[:19134])
f.close()

f=open("test_fold2.txt","w")
f.writelines(d1[19134:])
f.close()
