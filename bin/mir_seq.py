from __future__ import division
import shlex, subprocess
import commands
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

sp=raw_input("sp: ")
pre=raw_input("pre-file-name: ")

f= open("genomes/"+sp+"/mir/mir.gff3","r")
d1=f.readlines()
f.close()

d="\n".join(d1).split("#\n\n")[-1].split("\n\n")



record=[]
chrtmp=""
for x in d:
    t=x.split("\t")
    chr=t[0]
    start=int(t[3])
    end=int(t[4])
    st=t[6]
    name=t[-1].split("Name=")[-1].split(";")[0]
    if ("miRNA_primary_transcript" in x)==False:
        if chr!=chrtmp:
            f=open("genomes/"+sp+"/"+pre+chr+".fa")
            d1=SeqIO.read(f,"fasta")
            f.close()
            chrtmp=chr
        if st=="+":
            seq=d1.seq[start-1:end]
        else:
            seq=(d1.seq[start-1:end]).reverse_complement()
        content=SeqRecord(seq,name,'','')
        record.append(content)

f = open("genomes/"+sp+"/mir/miR.fasta", "w")
SeqIO.write(record, f, "fasta")
f.close()



