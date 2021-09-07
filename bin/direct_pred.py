#!/usr/bin/env python

from __future__ import division
import shlex, subprocess
import commands
import re
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

def disect(seq_name,seq_chr,chromosome0,seq_start,seq_end,seq_st,sequence,mir_length_container):        
    my_dna = Seq(sequence, generic_dna)
    sequence1=my_dna.transcribe()
    seq_frag=[]
    if mir_length_container+10<len(str(sequence1)):
        a=0
        while a<len(sequence1)-100:
            win=mir_length_container
            t=str(sequence1[a:a+win]).upper()
            name=">"+seq_name+";"+chromosome0+":"+str(seq_start+a)+"-"+str(seq_start+a+win-1)+":"+seq_st+"\n"
            seq_frag.append(name)
            seq_frag.append(t+"\n")
            a+=10
    else:
        name=">"+seq_name+";"+chromosome0+":"+str(seq_start)+"-"+str(seq_end)+":"+seq_st+"\n"
        seq_frag.append(name)
        seq_frag.append(str(sequence1).upper()+"\n")

    return seq_frag

def bp(struct):
    bp_list=[]

    pos=0
    left=[]
    right=[]
    for dot in struct:
        if dot=="(":
            left.append(pos)
        if dot==")":
            right.append(pos)
        pos+=1
    
    bp_p=0
    for left_pos in left:
        bp_list.append((left[bp_p],right[bp_p]))
        bp_p+=1

    return bp_list

def dinuc(sequence,dn_list):
    pos=0
    dn_c=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    while pos<len(sequence)-1:
        dn_c[dn_list.index(sequence[pos:pos+2])]+=1
        pos+=1
    dn_pos=1
    content=""
    for dn_val in dn_c[:16]:
        content=content+" "+str(dn_pos)+":"+str(dn_val/len(sequence))
        dn_pos+=1

    return content

def gu(sequence, bps_list):
    gu_c=0
    for pairs in bps_list:
        if sequence[pairs[0]]=="G" and sequence[pairs[1]]=="U":
            gu_c+=1
        if sequence[pairs[1]]=="G" and sequence[pairs[0]]=="U":
            gu_c+=1

    return gu_c

def triplet(sequence,struct,asb_list):
    mtx=[]
    m1x=[]
    i = 0
    while i<(len(sequence)-2):
        if struct[i]==")":
            tmp1="("
        else:
            tmp1=struct[i]
        if struct[i+1]==")":
            tmp2="("
        else:
            tmp2=struct[i+1]
        if struct[i+2]==")":
            tmp3="("
        else:
            tmp3=struct[i+2]
        ele=tmp1+tmp2+tmp3
        mtx.append(ele)
        m1x.append(sequence[i+1])
        i+=1

    triplet_fea=""
    fx = 21
    for feas in asb_list:
        i = 0
        val=0
        while i<len(mtx):
            if feas[1:] == mtx[i] and feas[0]==m1x[i]:
                val += 1
            i += 1
        fx += 1
        fval = val/len(sequence)
        triplet_fea=triplet_fea+" "+str(fx)+":"+str(fval)

    
    return triplet_fea
    

def features(sequence,struct,energy,dn_list,asb_list,group):  
    #find loops
    re1=re.compile("\(\.*\)")
    loop=re1.findall(struct)
    # number of loops
    n_loop=len(loop)
    # max loop size, feature 17
    max_loop=len(max(loop))-2
    # base pair list
    bps=bp(struct)
    # total base pairs
    tot_bp=len(bps)
    # avg_bp_stem, feature 18
    avg_bp_stem=tot_bp/len(sequence)
    # MFEI4, feature 19
    mfei4=float(energy)/tot_bp
    # dp/n_loop, feature 20
    fea20=tot_bp/len(sequence)/n_loop
    # GU property, feature 21
    fea21=gu(sequence, bps)/tot_bp

    # write all features
    f_content=group+dinuc(sequence,dn_list)+" 17:"+str(max_loop)+" 18:"+str(avg_bp_stem)+" 19:"+str(mfei4)+" 20:"+str(fea20)+" 21:"+str(fea21)+triplet(sequence,struct,asb_list)+" 54:"+energy+"\n"
    
    return f_content


root="/home/shared"

dn=["AA","AU","AC","AG","UA","UU","UC","UG","CA","CU","CC","CG","GA","GU","GC","GG","--","-A","-U","-G","-C","A-","U-","G-","C-"]
asb=["(((",".((","(.(","((.","..(",".(.","(..","..."]
nucl=['A','U','C','G']
fea=[]
for a in nucl:
    for b in asb:
        fea.append(a+b)

f=open("dilin_seq.txt","r")
d=f.readlines()
f.close()

seq=d[0]

dfrag=disect("dilin","chrN","chrN",0,4842,"+/-",seq,87)


f=open("dilin_testseq.txt","w")
f.writelines(dfrag)
f.close()

cmd="less "+root+"/mireval/uploads/dilin_testseq.txt|RNAfold --noPS > "+root+"/mireval/uploads/dilin_testfold.txt"
sta,out=commands.getstatusoutput(cmd)

f = open(root+"/mireval/uploads/dilin_testfold.txt","r")
data=f.readlines()
f.close()

f = open(root+"/mireval/uploads/dilin_test_features.txt",'w')
n=0
while n<len(data):
    seq4svm =data[n+1][:-1].upper()
    fold=data[n+2].split(" ")[0]
    mfe=data[n+2].split(" ")[-1].split(")")[0].split("(")[-1]
    head_mark="0"
    temp_svm=features(seq4svm,fold,mfe,dn,fea,head_mark)
    f.write(temp_svm)
    n+=3
f.close()

cmd="/home/shared/mireval/bin/svm_light/svm_classify "+root+"/mireval/uploads/dilin_test_features.txt /home/shared/mireval/ref/svm/model.txt "+root+"/mireval/uploads/dilin_res_svm.txt"
sta,out=commands.getstatusoutput(cmd)
f=open(root+"/mireval/uploads/dilin_res_svm.txt","r")
svm_res=f.readlines()
f.close()


pos_res=[]
i=0
for x in svm_res:
    if float(x[:-1])>0:
        pos_res.append(data[i])
        pos_res.append(data[i+1])
        pos_res.append(data[i+2])
    i+=3

f=open("pos_res","w")
f.writelines(pos_res)
f.close()
