from __future__ import division
import shlex, subprocess
import commands
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.Blast import NCBIXML



def conserve_score(root_name,chromosome,start_cor,end_cor,species,method,group):
    run="less "+root_name+"/web/genomes/"+species+method+"/chr"+chromosome+".phastCons46way."+group+".wigFix|grep fix > temp.txt "
    sta,out=commands.getstatusoutput(run)   
    
    files=open(root+"/web/genomes/"+species+method+"temp.txt","r")
    step_list=files.readlines()
    files.close()
    tracer=""
    for steps in step_list:
        step_start=int(steps.split("start=")[1].split(" step")[0])
        if start_cor >= step_start:
            tracer = steps
        else:
            break


    files=open(root_name+"/web/genomes/"+species+method+"/chr"+chromosome+".phastCons46way."+group+".wigFix","r")
    score_list=files.readlines()
    files.close()
    
    if tracer == "":
        i=0
    else:    
        i=score_list.index(tracer)
    tmp=start_cor
    h=0
    plot=[]
    
    while i < len(score_list):
        if "fix" in score_list[i]:
            point=int(score_list[i].split("start=")[1].split(" step")[0])
            if h == 1 and point <= end_cor:
                plot=plot+["0"]*(point-tmp)
                tmp=point
            if point> end_cor:
                plot=plot+["0"]*(end_cor-tmp+1)
                break
            
        else:
            if point == tmp:
                plot.append(score_list[i][:-1])
                tmp += 1
                h=1
            if point == end_cor:
                break
            point+=1 
         
        i+=1    
        
    run="rm temp.txt"
    sta,out=commands.getstatusoutput(run)     
    return plot



def neighbour(root_name,chromosome,start_cor,end_cor,strand,species):
    files=open(root_name+"/web/genomes/"+species+"/mir/mir.gff3","r")
    mir_list1=files.readlines()
    files.close()
    mir_list="\n".join(mir_list1).split("#\n\n")[-1].split("\n\n")

    upstream=[]
    for premir in mir_list:
        if ("miRNA_primary_transcript" in premir)==False:
            mirchr="chr"+premir.split("\t")[0]
            mirstart=int(premir.split("\t")[3])
            mirend=int(premir.split("\t")[4])
            mirname=premir.split("\t")[-1].split("Name=")[-1].split(";")[0]
            mirst=premir.split("\t")[6]
            if mirchr==chromosome:
                if (start_cor - 5000) < mirend and mirend < star_cor:
                    len1=mirend-mirstart+1
                    dis=star_cor-mirend
                    upstream.append(mirname+"_"+mirst+"_"+str(dis)+"_"+len1+"_upstream")
                if (end_cor + 5000) > mirstart and mirstart > end_cor:
                    len1=mirend-mirstart+1
                    dis=mirstart-end_cor
                    upstream.append(mirname+"_"+mirst+"_"+str(dis)+"_"+len1+"_downstream")
                if mirend < end_cor and mirstart > start_cor:
                    len1=mirend-mirstart+1
                    dis=mirstart-start_cor
                    upstream.append(mirname+"_"+mirst+"_"+str(dis)+"_"+len1+"_within")
    return upstream


def disect(root_name,seq_name,seq_chr,q_seq_start,q_seq_end,seq_st,sequence,mir_length_container):        
    seq_frag=[]
    if max(mir_length_container)+10<len(sequence):
        a=0
        while a<len(sequence)-100:
            win=random.sample(mir_length_container,1)[0]
            t=str(sequence[a:a+win])
            name=">"+seq_name+":"+seq_chr+":"+str(q_seq_start+a)+"-"+str(q_seq_start+a+win-1)+":"+seq_st
            seq_frag.append(t)
            a+=10
    else:
        seq_frag.append(">"+seq_name+":"+seq_chr+":"+str(q_seq_start)+"-"+str(q_seq_end)+":"+seq_st+"\n")
        seq_frag.append(sequence+"\n")

    return seq_frag 


def polii(seq_name,chromosome,start_cor,end_cor,st,sequence):
    seq_polii=[]
    seq_polii.append(">"+seq_name+":"+chromosome+":"+str(start_cor)+"-"+str(end_cor)+":"+st)
    seq_polii.append(sequence)
    
    return seq_polii



#<Main>

root="/home/dadi/projects/web/"
f=open(root+"/web/res/record_log.txt","r")
record_step=f.readlines()
f.close()




res_no=str(int(record_step[0][:-1])+1)
cmd="mkdir "+root+"/web/res/"+res_no
sta,out=commands.getstatusoutput(cmd)


seq_file=raw_input("Input File: ")

sp=raw_input("Species Name: ")

classes=raw_input("Class: ")    # primates, worms, vertebrate, dog, horse, prozeon

hit_thresold=raw_input("Hit Thresold (%): ")
hit_thres=float(hit_thresold)/100





# command to run blast

cmd="blastall -p blastn -d "+root+"/web/genome/"+species+"/"+database+"_ref.fa -i "+root+"/web/res/"+seq_file+"-m 7 -o "+root+"/web/res/"+res_no+"/"+"temp.txt"
sta,out=commands.getstatusoutput(cmd)


# get best-hit of blast and its coordinates for conservation studies


blast_records = NCBIXML.parse(open(root+"/web/res/"+res_no+"/temp.txt"))
#open("res/"+res_no+"/temp.txt").close()

dblast=[]
dfrag=[]
dpolii=[]
blast_no=0
for blast_record in blast_records:
    name=blast_record.query
    len=blast_record.query_letters
    dblast.append(name+" "+"begin\n\n\n")

    f=open(root+"/web/genomes/"+sp+"/mir/mir.gff3","r")
    mir_info=f.readlines()
    f.close()
    
    mir_length=[]
    for mir in mir_info[:-1]:
        if ("primary" in mir) == True and mir[0]!="#":
            t=mir.split("\t")
            mir_length.append(int(t[4])-int(t[3])+1)  


    
    for alignment in blast_record.alignments:
        hsp=alignment.hsp[0]
        if hsp.identities<len*hit_thres:
            dblast.append("Blast: "+"Not mappable\n")
            break
        else:
            blast_no+=1
            blast_name=name+"_"+str(blast_no)
            query_start=hsp.query_start
            query_end=hsp.query_end
            start=hsp.sbjct_start
            end=hsp.sbjct_end
            if hsp.frame==(1,1):
                st="+"
            else:
                st="-"
            chr=(alignment.title).split("chromosome ")[1].split(",")[0]
            seq=str(hsp.sbjct).upper()
            content=blast_name+":"+"chr"+chr+":"+str(start)+"-"+str(end)+":"+st+"\n"                          # blast res                                                              
            plot1=conserve_score(root,chr,start,end,sp,"/cons","")                                            # conservation res
            plot2=conserve_score(root,chr,start,end,sp,"/shadow",classes)                                     # shadow res
            nearby=neighbour(root,chr,start,end,st,sp)
            dfrag=dfrag+disect(root,blast_name,chr,query_start,query_end,st,seq,mir_length)
            dpolii=dpolii+polii(root,blast_name,start,end,st,seq)
                   

    dblast.append(name+" "+"end\n\n\n")





f=open(root+"/web/res/"+res_no+"/res_blast.txt","w")
f.writelines(dblast)
f.close()

f=open(root+"/web/HeteroMirPred/res_frag_"+res_no+".fasta","w")
f.writelines(dfrag)
f.close()




# Find + prediction by HeteroMirPred

cmd="perl HeteroMirPred.pl res_frag_"+res_no+".fasta"
sta,out=commands.getstatusoutput(cmd)

f=open(root+"/web/HeteroMirPred/res_frag_"+res_no+".fasta","r")
d=f.readlines()
f.close()

preds=out.split("##\n")[-1].split("\n")


positive=[]
pred_no=1
fasta_no=0
for pred in preds:
    if int(pred.split(" ")[-1])>=0.5:
        positive.append((fasta_no,pred_no))
        pred_no+=1
    fasta_no+=1


struct2nd_rec=out.split("zdata\n")[1:]


pred_res=[]

for x in positive:
    pred_res.append(d(x[0]*2)[:-1]+":"+str(x[1])+"\n")
    pred_res.append(struct2nd_rec(x[0]*2).split("\n")[0]+"\n")
    pred_res.append(struct2nd_rec(x[0]*2).split("\n")[1]+"\n")
    
f=open(root+"/web/res/"+res_no+"/res_pred.txt","w")
f.writelines(pred_res)
f.close()

cmd="rm "+root+"/web/HeteroMirPred/res_frag_"+res_no+".fasta"
sta,out=commands.getstatusoutput(cmd)


# Promoter Finder
f=open(root+"/web/res_4fimo_"+res_no+".txt","w")
f.writelines(dpolii)
f.close()

cmd="fimo "+root+"/web/pol2.meme "+root+"/web/res_4fimo_"+res_no+".txt"
sta,out=commands.getstatusoutput(cmd)

cmd="mv "+root+"/web/fimo_out/fimo.txt "+root+"/web/res/"+res_no+"/res_fimo.txt"
sta,out=commands.getstatusoutput(cmd)

cmd="rm "+root+"/web/res_4fimo_*"
sta,out=commands.getstatusoutput(cmd)

cmd="rm "+root+"/web/fimo_out/*"
sta,out=commands.getstatusoutput(cmd)



# merge result and create a genbank file



    


f=open(root+"/web/record_log.txt","w")
f.writelines([res_no+"\n"])
f.close()
