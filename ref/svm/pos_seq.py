#!/usr/bin/python


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

root="/home/dadi/projects"

sp_container=["homo_sapiens", "mus_musculus", "bos_taurus", "caenorhabditis_briggsae", "caenorhabditis_elegans", "canis_familiaris", "drosophila_melanogaster", "danio_rerio", "gallus_gallus", "xenopus_tropicalis", "rattus_norvegicus", "pan_troglodytes", "macaca_mulatta", "monodelphis_domestica", "fugu_rubripes"]

d1=[]
d2=[]
for sp in sp_container:
    print sp
    new1=[]
    new2=[]
    mf=open(root+"/web/genomes/"+sp+"/mir/mirseq.txt","r")
    d=mf.readlines()
    mf.close()
    times=int(int(len(d)/2)*0.9)
    list0=random.sample(range(0,int(len(d)/2)),int(len(d)/2))
    list1=list0[:times]
    list2=list0[times:]
    for y in list1:
        if d[2*y][0]==">" and d[2*y+1][0] != ">" and d[2*y+1][0] != "E":
            new1.append(d[2*y])
            new1.append(d[2*y+1])
    for z in list2:
        if d[2*z][0]==">" and d[2*z+1][0] != ">" and d[2*z+1][0] != "E":
            new2.append(d[2*z])
            new2.append(d[2*z+1])
    d1=d1+new1
    d2=d2+new2

print len(d1)
print len(d2)

f=open("train_pos_seq_rep1_all.txt","w")
f.writelines(d1)
f.close()

f=open("test_pos_seq_rep1_all.txt","w")
f.writelines(d2)
f.close()
