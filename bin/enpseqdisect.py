from __future__ import division
import linecache
import shlex, subprocess
import commands
import re
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


f=open("bonafold.txt","r")
da=f.readlines()
f.close()

l=[]
i=0
while i<len(da):
    l.append(len(da[i+1]))
    i+=3

f=open("rn45s.fastq","r")
data=f.readlines()
f.close()


d=[]
i=0
g=0
while i<len(data):
    g+=1
    seq=data[i+1]
    a=0
    while a<len(seq)-100:
        win=random.sample(l,1)[0]
        t=str(seq[a:a+win])
        name=data[i][:-2]+" "+str(a)+"-"+str(a+win-1)+"\n"
        d.append(name)
        d.append(t+"\n")
        a+=10

    i+=2



f=open("rn45sdisect.txt","w")
f.writelines(d)
f.close()

