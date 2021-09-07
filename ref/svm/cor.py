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

root="/home/shared"

sp_container=["homo_sapiens", "mus_musculus", "bos_taurus", "caenorhabditis_briggsae", "caenorhabditis_elegans", "canis_familiaris", "drosophila_melanogaster", "danio_rerio", "gallus_gallus", "xenopus_tropicalis", "rattus_norvegicus", "pan_troglodytes", "macaca_mulatta", "monodelphis_domestica", "fugu_rubripes"]


for sp in sp_container:
    print sp
    new=[]
    chrf=open(root+"/mireval/genomes/"+sp+"/chrom.sizes","r")
    d0=chrf.readlines()
    chrf.close()
    mf=open(root+"/mireval/genomes/"+sp+"/mir/mir.gff3","r")
    d=mf.readlines()
    mf.close()
    for mir in d:
        if ("miRNA_primary_transcript" in mir) and mir[0]!="#":
            t=mir.split("\t")
            if t[0][:3]!="chr":
                mirchr="chr"+ t[0]
            else:
                mirchr=t[0]
            if mirchr[-2:]==".1":
                mirchr=mirchr[:-2]
            mirstart=t[3]
            mirend=t[4]
            if t[6]=="-":
                mirst="minus"
            else:
                mirst="plus"
            mirname=t[-1].split("Name=")[-1].split(";")[0][:-1]
            #print mirchr
            for chrname in d0:
                if mirchr==chrname.split("\t")[0] or mirchr[3:]==chrname.split("\t")[0]:
                    chr0=chrname.split("\t")[0]
                    break
            cmd="blastdbcmd -db "+root+"/mireval/genomes/"+sp+"/blastdb -dbtype nucl -entry "+chr0+" -strand "+mirst+" -range "+mirstart+"-"+mirend+" -outfmt %s"
            #print cmd
            sta,out=commands.getstatusoutput(cmd)
            new.append(">"+mirname+"\n")
            new.append(out+"\n")
    f=open(root+"/mireval/genomes/"+sp+"/mir/mirseq.txt","w")
    f.writelines(new)
    f.close()

