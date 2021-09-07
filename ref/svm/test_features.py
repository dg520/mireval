#!/usr/bin/python

from __future__ import division
import linecache
import shlex, subprocess
import commands
import re


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
    dn_c=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
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
    fx = 24
    #fx=21
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
    
def bulge(struct,loop_no):

    re3=re.compile("\(\.\.*\(")
    buldge_l=re3.findall(struct)
    if buldge_l==[]:
        max_buldge_l=0
    else:
        max_buldge_l=len(max(buldge_l))-2
    
    re4=re.compile("\)\.\.*\)")
    buldge_r=re4.findall(struct)
    if buldge_r==[]:
        max_buldge_r=0
        tot_br=0
    else:
        max_buldge_r=len(max(buldge_r))-2

    if loop_no==0:
        max_buldge=len(struct)

    else:
        if loop_no==1 or loop_no==2:
            max_buldge=max_buldge_r+max_buldge_l
            #max_buldge=max(max_buldge_r,max_buldge_l)
        else:
            max_buldge=loop_no*max(max_buldge_r,max_buldge_l)


    return max_buldge


def features(sequence,struct,energy,dn_list,asb_list,group):  
    #find loops
    re1=re.compile("\(\.*\)")
    loop=re1.findall(struct)

    re2=re.compile("\.*")
    buldge=re2.findall(struct)
    # total buldges, feature 24
    tot_buldge=len("".join(buldge))
    # maximum buldge size, feature 23
    max_bulge=bulge(struct,len(loop))


    # base pair list
    bps=bp(struct)
    
    # total base pairs
    if len(loop)==0 or len(loop)==1:
        tot_bp=len(bps)
    else:
        re5=re.compile("\)\.*\(")
        inner_bp=re5.findall(struct)
        max_inner_bp=max(inner_bp)
        brak_r=0
        brak_l=0
        for brak in max_inner_bp:
            if brak=="(":
                brak_l+=1
            if brak==")":
                brak_r+=1
        tot_bp=max(brak_r,brak_l)
    # number of loops, feature 22 and max loop size, feature 17 and dp/n_loop, feature 20
    if loop != []:
        loop_size=0
        n_loop=len(loop)
        if n_loop==1 or n_loop==2:
            n_loop_res=1
        else:
            n_loop_res=0
        max_loop=len(max(loop))-2
        fea20=tot_bp/len(sequence)/n_loop
    else:
        n_loop=0
        n_loop_res=0
        max_loop=0
        fea20=0

    
    
    # avg_bp_stem, feature 18
    avg_bp_stem=tot_bp/len(sequence)
    # MFEI4, feature 19 and GU property, feature 21
    if tot_bp !=0:
        mfei4=float(energy)/tot_bp
        fea21=gu(sequence, bps)/tot_bp
    else:
        mfei4=0
        fea21=0

    # write all features
    f_content=group+dinuc(sequence,dn_list)+" 17:"+str(max_loop)+" 18:"+str(avg_bp_stem)+" 19:"+str(mfei4)+" 20:"+str(fea20)+" 21:"+str(fea21)+" 22:"+str(n_loop_res)+" 23:"+str(max_bulge)+" 24:"+str(tot_buldge)+triplet(sequence,struct,asb_list)+" 57:"+energy+"\n"
    #f_content=group+dinuc(sequence,dn_list)+" 17:"+str(max_loop)+" 18:"+str(avg_bp_stem)+" 19:"+str(mfei4)+" 20:"+str(fea20)+" 21:"+str(fea21)+triplet(sequence,struct,asb_list)+" 54:"+energy+"\n"
    
    return f_content
    
    
    
    


dn=["AA","AU","AC","AG","UA","UU","UC","UG","CA","CU","CC","CG","GA","GU","GC","GG"]
asb=["(((",".((","(.(","((.","..(",".(.","(..","..."]
nc=['A','U','C','G']
fea=[]
for a in nc:
    for b in asb:
        fea.append(a+b)


f = open("test_fold_rep1.txt","r")
data=f.readlines()
f.close()




f = open("test_features_rep1_822.txt",'w')
n=0
while n<len(data):
    seq =data[n+1][:-1].upper()
    fold=data[n+2].split(" ")[0]
    mfe=data[n+2].split(" ")[-1].split(")")[0].split("(")[-1]


    head="0"

    temp=features(seq,fold,mfe,dn,fea,head)
    f.write(temp)
    n+=3
    #print n/3
f.close()