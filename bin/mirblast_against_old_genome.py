from __future__ import division
import shlex, subprocess
import commands
from Bio.Blast import NCBIXML



root="/home/shared/mireval/genomes/"



for sp in ["bos_taurus","caenorhabditis_briggsae","caenorhabditis_elegans"]:
    # seq in new genome
   # print sp
   # print "fetching seq in new genome..."
   # f=open(root+sp+"/mir/mir.gff3","r")
   # mir_info=f.readlines()
   # f.close()
   # premir=[]
   # for mir in mir_info:
   #     if ("primary" in mir) == True and mir[0]!="#":
   #         t=mir.split("\t")
   #         if t[6]=="+":
   #             st="plus"
   #         else:
   #             st="minus"
   #         cmd="/usr/bin/blastdbcmd -db "+root+sp+"/current_genome/blastdb -dbtype nucl -entry "+t[0]+" -strand "+st+" -range "+t[3]+"-"+t[4]+" -outfmt %s"
   #         sta,out=commands.getstatusoutput(cmd)
   #         name=t[-1].split("=")[-1]
   #         premir.append(">"+name)
   #         premir.append(out+"\n")
   # f=open(root+sp+"/mir/primirseq.txt","w")
   # f.writelines(premir)
   # f.close()
    # blast against old genome
   # print "blasting against old genome..."
   # cmd="/usr/bin/blastn -db /home/shared/mireval/genomes/"+sp+"/blastdb -query /home/shared/mireval/genomes/"+sp+"/mir/primirseq.txt -outfmt 5 -out /home/shared/mireval/genomes/"+sp+"/mir/blastres.txt"
   # sta,out=commands.getstatusoutput(cmd)
    # pick out hits
   # print "picking out best hits..."
   # blast_records = NCBIXML.parse(open(root+sp+"/mir/blastres.txt"))
   # d2=[]
   # for blast_record in blast_records:
   #     blast_no=0
   #     name_hit=str(blast_record.query)
   #     hit_len=blast_record.query_letters
   #     hn=0
   #     for alignment in blast_record.alignments:
   #         hsp=alignment.hsps[0]
   #         if hsp.identities+1==hit_len or hsp.identities==hit_len:
   #             hn+=1
   #             chr=str((alignment.title).split("chromosome ")[-1].split(",")[0]).split(" No definition line")[0]
   #             if hsp.frame==(1,1):
   #                 st_hit="+"
   #             else:
   #                 st_hit="-"
   #             d2.append(chr+"\t"+str(min(int(hsp.sbjct_start),int(hsp.sbjct_end)))+"\t"+str(max(int(hsp.sbjct_start),int(hsp.sbjct_end)))+"\t"+name_hit+"\t"+st_hit+"\n")
   #     if hn>=2 or hn==0:
   #         print name_hit, hn
         
   # f=open(root+sp+"/mir/primir_old_cor.txt","w")
   # f.writelines(d2)
   # f.close()

    # prepare new gff3
    print "preparing new gff3"
    f=open(root+sp+"/mir/mir.gff3","r")
    mir_info=f.readlines()
    f.close()
    d3=[]
    head=[]
    content=[]
    sub=[]
    for x in mir_info:
        if x[0]=="#":
            head.append(x)
        elif ("primary" in x):
            content.append(sub)
            sub=[]
            sub.append(x)
        else:
            sub.append(x)
    content.append(sub)
    content=content[1:]
    
 
    f=open(root+sp+"/mir/primir_old_cor.txt","r")
    oc=f.readlines()
    f.close()

    if len(content)!=len(oc):
        print "Error!"
        break
    new=[]
    i=0
    while i<len(content):
        t1=content[i][0].split("\t")
        t2=content[i][0].split("\t")
        t3=oc[i].split("\t")
        start=int(t3[1])
        end=int(t3[2])
        t2[0]=t3[0]
        t2[3]=t3[1]
        t2[4]=t3[2]
        t2[6]=oc[i][-2]
        new.append("\t".join(t2))
        for x in content[i][1:]:
            tt=x.split("\t")
            diff1=int(tt[3])-int(t1[3])
            diff2=int(t1[4])-int(tt[4])
            tt[0]=t3[0]
            tt[3]=str(start+diff1)
            tt[4]=str(end-diff2)
            tt[6]=oc[i][-2]
            new.append("\t".join(tt))
        i+=1
    d3=head+new
    
    f=open(root+sp+"/mir/mir_new.txt","w")
    f.writelines(d3)
    f.close()
