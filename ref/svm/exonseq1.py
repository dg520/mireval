from __future__ import division
import linecache
import shlex, subprocess
import commands
import re
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


root="/home/dadi/projects"
sp="mus_musculus"

l=[83,84]


f=open("xaa","r")
data=f.readlines()
f.close()

data=data[1:]

d=[]
c=0
no=""
g=0
theno=1
for x in data:
    temp=x.split("\t")
    if "_" in temp[0]:
        break
    else:
        i=0
        for y in temp[2].split(",")[:-1]:
            g=g+1
            start=int(y)
            end=int(temp[3].split(",")[:-1][i])

            if (end-start+1)>=84:
                if temp[1]=="+":
                    st="plus"
                else:
                    st="minus"
                cmd="blastdbcmd -db "+root+"/web/genomes/"+sp+"/blastdb -dbtype nucl -entry "+temp[0]+" -strand "+st+" -range "+str(start)+"-"+str(end)+" -outfmt %s"
                sta,out=commands.getstatusoutput(cmd)
                seq=out
                #print out
                a=0
                while a<len(seq)-94:
                    win=random.sample(l,1)[0]
                    t=str(seq[a:a+win])
                    name=">"+temp[0]+"\t"+str(start+a)+"\t"+str(start+a+win-1)+"\t"+temp[1]+"\t"+str(g)+"\n"
                    d.append(name)
                    d.append(t+"\n")
                    a+=10
            i+=1
        c+=1
    
        if c>147:
            print 1*theno,"%"
            #print d[-2:]
            c=0
            theno+=1

f=open("exonseq83random1.txt","w")
f.writelines(d)
f.close()


