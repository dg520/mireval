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
    
    files=open("temp.txt","r")
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


def disect(root_name,sequence,species):
    f=open(root_name+"/web/genomes/"+species+"/mir/mir.gff3","r")
    d=f.readlines()
    f.close()
    
    mir_length=[]
    for mir in d[:-1]:
        if "primary" in mir == True:
            t=mir.split("\t")
            mir_length.append(int(t[4])-int(t[3])+1)
   
        
    seq_frag=[]    
    a=0
    while a<len(sequence)-100:
        win=random.sample(mir_length,1)[0]
        t=str(sequence[a:a+win])
        seq_frag.append(t)
        a+=10
    
    res_list=[max(mir_length)]+seq_frag=[]
    return res_list 



    
#def mirblast(root_name,seq_file,species,database):
#    run="blastall -p blastn -d "+root_name+"/web/genome/"+species+"/mir/miR.fasta -i "+root_name+"/web/genome/"+species+"/"+seq_file+" -o temp.txt"
#    sta,out=commands.getstatusoutput(run)


#    files=open("temp.txt","r")
#    blast_list=files.readlines()
#    files.close()    
    

#    dblast=[]
#    
#    for record in blast_list:
#        length=record.split("\n")[1].split("(")[-1].split(" letters")[0]
#        name=record.split("\n")[0]

#        dtem=record.split(">")[1:]
#        dblast.append(name+" "+"begin\n\n\n")

#        for lines in dtem:
#            id= lines.split("Identities = ")[-1].split("/")[0]
#            if int(id) < length*99%:
#                dblast.append("Blast: "+"Not mappable\n")
#                break
#            else:
#                n1=lines.split("Sbjct: ")[1].split(" ")[0]
#                n2=lines.split(" ")[-1].split("\n\n\n")[0]
#                start=min(int(n1),int(n2))
#                end=max(int(n1),int(n2))
#                chr=lines.split("chromosome ")[1].split(",")[0]
#                st=lines.split("Plus / ")[1].split("\n")[0]
#                content="Blast: "+"chr"+chr+":"+str(start)+"-"+str(end)+":"+st+"\n"                     # blast res
#                plot1=conserve_score(chr,start,end,sp,"/cons","")
#                plot2=conserve_score(chr,start,end,sp,"/shadow",classes) 
#    dblast.append(name+" "+"end\n\n\n")
#    return dblast





#<Main>
f=open(root_name+"/web/res/record_log.txt","r")
record_step=f.readlines()
f.close()




res_no=str(int(record_step[0][:-1])+1)
cmd="mkdir res/"+res_no
sta,out=commands.getstatusoutput(cmd)


seq_file=raw_input("Input File: ")
root="/home/dadi/projects/web/"
sp=raw_input("Species Name: ")

classes=raw_input("Class: ")    # primates, worms, vertebrate, dog, horse, prozeon

hit_thresold=raw_input("Hit Thresold (%): ")
hit_thres=float(hit_thresold)/100





# command to run blast

cmd="blastall -p blastn -d "+root_name+"/web/genome/"+species+"/"+database+"_ref.fa -i "+root_name+"/web/res/"+seq_file+"-m 7 -o "+root_name+"/web/res/"+res_no+"/"+"temp.txt"
sta,out=commands.getstatusoutput(cmd)


# get best-hit of blast and its coordinates for conservation studies


blast_records = NCBIXML.parse(open("res/"+res_no+"/temp.txt"))
#open("res/"+res_no+"/temp.txt").close()

dblast=[]
for blast_record in blast_records:
    name=blast_record.query
    len=blast_record.query_letters
    dblast.append(name+" "+"begin\n\n\n")

    for alignment in blast_record.alignments:
        hsp=alignment.hsp[0]
        if hsp.identities<len*hit_thres:
            dblast.append("Blast: "+"Not mappable\n")
            break
        else:
            start=hsp.query_start
            end=hsp.query_end
            if hsp.frame==(1,1):
                st="+"
            else:
                st="-"
            chr=(alignment.title).split("chromosome ")[1].split(",")[0]
            seq=str(hsp.sbjct).lower()
            content="Blast: "+"chr"+chr+":"+str(start)+"-"+str(end)+":"+st+"\n"                          # blast res                                                              
            plot1=conserve_score(root,chr,start,end,sp,"/cons","")                                       # conservation res
            plot2=conserve_score(root,chr,start,end,sp,"/shadow",classes)                                # shadow res
            nearby=neighbour(root,chr,start,end,st,sp)
            seq_disect_temp=disect(seq)
            mir_max=seq_disect_temp[0]
            if len(seq)<=seq_disect_temp[0]+10:
                #mirpred directly
            else:
                for x in seq_disect_temp[1:]:
                    #mirpred directly

    dblast.append(name+" "+"end\n\n\n")











#f=open("res/"+res_no+"/"+"temp.txt","r")                   # File blast name
#d=f.readlines()
#f.close()


#dblast=[]


#d1="".join(d).split("  Database:")[0].split("Query= ")[1:]



#for x in d1:
#    len=x.split("\n")[1].split("(")[-1].split(" letters")[0]
#    name=x.split("\n")[0]

#    dtem=x.split(">")[1:]
#    dblast.append(name+" "+"begin\n\n\n")

#    for y in dtem:
#        id= y.split("Identities = ")[-1].split("/")[0]
#        if int(id) < len*99%:
#            dblast.append("Blast: "+"Not mappable\n")
#            break
#        else:
#            n1=y.split("Sbjct: ")[1].split(" ")[0]
#            n2=y.split(" ")[-1].split("\n\n\n")[0]
#            start=min(int(n1),int(n2))
#            end=max(int(n1),int(n2))
#            chr=y.split("chromosome ")[1].split(",")[0]
#            st=y.split("Plus / ")[1].split("\n")[0]
#
#            seq=

#            content="Blast: "+"chr"+chr+":"+str(start)+"-"+str(end)+":"+st+"\n"                          # blast res                                                              
#            plot1=conserve_score(root,chr,start,end,sp,"/cons","")                                       # conservation res
#            plot2=conserve_score(root,chr,start,end,sp,"/shadow",classes)                                # shadow res
#            nearby=neighbour(root,chr,start,end,st,sp)
                
#    dblast.append(name+" "+"end\n\n\n")


f=open(root_name+"/web/res/"+res_no+"/res.txt","w")
f.writelines(dblast)
f.close()


# create a genbank file



    
    


f=open("record_log.txt","w")
f.writelines([res_no+"\n"])
f.close()
