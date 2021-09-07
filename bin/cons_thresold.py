from __future__ import division
import os, sys, getopt
import shlex, subprocess
import commands
import numpy
import math


def usage():
    print """
USAGE: ./cons_thresold.py -v <chrom_size_file> -C <conservertion_score.bw> -S <shadowing_score.bw> -M <mir.gff3>

Options:
-v              Chromosome-size file. [Compulsory]
-C              Bigwig file that contains multiple aligment scores (please read TUTOR). [Compulsory]
-S              Bigwig file that contains multiple aligment scores for closely related species (please read TUTOR). [Optional]
-M              GFF3 file annotating miRNA in miRBase format. [Compulsory]
-h              Print this help
--help

"""
    return 1

def par(argv):
    #animal=""
    cz_loc=""
    cons_loc=""
    shad_loc=""
    mir_loc=""
    try:
        opts, args = getopt.getopt(argv,"v:C:S:M:",["vfile=","Cfile=","Sfile=","Mfile="])
        #print opts, args
    except getopt.GetoptError:
        usage() # -X : Promyelocyte_dominent_ir_boundary; -Y : Granulocyte_dominent_ir_boundary
        sys.exit(2)
    for opt, arg in opts:
        #print opt,arg
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-v", "--vfile"):
            cz_loc = arg
        elif opt in ("-C", "--Cfile"):
            cons_loc= arg
        elif opt in ("-S", "--Sfile"):
            shad_loc= arg
        elif opt in ("-M", "--Mfile"):
            mir_loc= arg
        
        else:
            usage()
            sys.exit(2)

    return (cz_loc, cons_loc, shad_loc, mir_loc)





def percentile(N, P):
    n = int(round(P * len(N) + 0.5))
    return N[n-1]




def cons_list(fn):
    f=open(fn,"r")
    score_list=f.readlines()
    f.close()
    ini=[]
    if score_list !=[]:
        for y in score_list:
            t2=y.split("\t")
            for x in range(0,(int(t2[2])-int(t2[1]))):
                ini.append(float(t2[-1][:-1]))
    
    return ini 
    


def thres_by_mir_mature(root_name,dir_mir,dir_cz,dir_cons,method):
    f=open(dir_mir,"r")
    mir_info=f.readlines()
    f.close()
    mir_length=0
    val=[]
    av=[]
    for mir in mir_info:
        if ("primary" in mir) == False and mir[0]!="#":
            t=mir.split("\t")
            temp_length=int(t[4])-int(t[3])+1
            mir_length=mir_length+temp_length
            f=open(dir_cz,"r")
            chr_info=f.readlines()
            f.close()
            chr=0
            for x in chr_info:
                t1=x.split("\t")
                if t[0].split(".")[0] == t1[0]:
                    chr=t1[0]
                    break
            if chr==0:
                for x in chr_info:
                    t1=x.split("\t")
                    if (t[0].split(".")[0] in t1[0]):
                        chr=t1[0]
                        break
            if chr!=0:
                #print mir
                cmd=root_name+"/mireval/bin/bigWigToBedGraph "+dir_cons+" -chrom="+chr+" -start="+t[3]+" -end="+t[4]+" score.txt"
                sta,out=commands.getstatusoutput(cmd)
                file_name="score.txt"
                temp_list=cons_list(file_name)
                #if len(temp_list)<temp_length:
                #    for dist in range(0,(temp_length-len(temp_list))):
                #        temp_list.append(0.0)
                val=val+temp_list
                if temp_list !=[]:
                    #av.append(str(numpy.mean(temp_list)))
                    if numpy.mean(temp_list)>0.1:
                        av.append(numpy.mean(temp_list))
    #val=sorted(val,key=float)
    #cons_25=percentile(val,0.25)
    #cons_75=percentile(val,0.75)

    #f=open("/home/shared/mireval/genomes/"+species+"/mir/4hist","w")
    #for x in val:
    #    f.write(str(x)+"\n")
    #f.close()
    cons_mean=numpy.mean(av)
    cons_std=numpy.std(val)
    #print method[:-1], cons_mean, cons_std, cons_25, cons_75
    print method, cons_mean

    run="rm score.txt"
    sta,out=commands.getstatusoutput(run)    

    #return 1
    return av



def thres_by_mir(root_name,dir_mir,dir_cz, dir_cons,method):
    f=open(dir_mir,"r")
    mir_info=f.readlines()
    f.close()
    mir_length=0
    val=[]
    av=[]
    for mir in mir_info:
        if ("primary" in mir) == True and mir[0]!="#":
            t=mir.split("\t")
            temp_length=int(t[4])-int(t[3])+1
            mir_length=mir_length+temp_length
            f=open(dir_cz,"r")
            chr_info=f.readlines()
            f.close()
            chr=0
            for x in chr_info:
                t1=x.split("\t")
                if t[0].split(".")[0] == t1[0]:
                    chr=t1[0]
                    break
            if chr==0:
                for x in chr_info:
                    t1=x.split("\t")
                    if (t[0].split(".")[0] in t1[0]):
                        chr=t1[0]
                        break
            if chr!=0:
                #print mir
                cmd=root_name+"/mireval/bin/bigWigToBedGraph "+dir_cons+" -chrom="+chr+" -start="+t[3]+" -end="+t[4]+" score.txt"
                sta,out=commands.getstatusoutput(cmd)
                file_name="score.txt"
                temp_list=cons_list(file_name)
                if len(temp_list)<temp_length:
                    for dist in range(0,(temp_length-len(temp_list))):
                        temp_list.append(0.0)
                val=val+temp_list
                if temp_list !=[]:
                    #av.append(str(numpy.mean(temp_list)))
                    if numpy.mean(temp_list)>0.1:
                        av.append(numpy.mean(temp_list))
                #else:
                    #av.append(0.0)
    #val=sorted(val,key=float)
    #cons_25=percentile(val,0.25)
    #cons_75=percentile(val,0.75)

    #f=open("/home/shared/mireval/genomes/"+species+"/mir/4hist","w")
    #for x in val:
    #    f.write(str(x)+"\n")
    #f.close()
    cons_mean=numpy.mean(av)
    cons_std=numpy.std(val)
    #print method[:-1], cons_mean, cons_std, cons_25, cons_75
    print method, cons_mean
    
    run="rm score.txt"
    sta,out=commands.getstatusoutput(run)

    #return 1
    return av


cwd = os.getcwd()
root="/".join(cwd.split("/")[:-1])

cz_dir, cons_dir, shad_dir, mir_dir =par(sys.argv[1:])


if cz_dir!="" and cons_dir!="" and mir_dir!="":
    t_av1=thres_by_mir_mature(root,mir_dir,cz_dir,cons_dir,"conservertion upper thresold: ")
    t_av3=thres_by_mir(root,mir_dir,cz_dir,cons_dir, "conservertion lower thresold: ")

        
    # shad threshold
    if shad_dir!="":
        t_av2=thres_by_mir_mature(root,mir_dir,cz_dir,shad_dir, "shadowing upper thresold: ")
        t_av4=thres_by_mir(root,mir_dir,cz_dir,shad_dir,  "shadowing lower thresold: ")
    else:
        print "shadowing upper thresold: ", "NA"
        print "shadowing lower thresold: ", "NA"


else:
    usage()
    if cz_dir=="":
        print "Error: Need a chromosome-size file."
    if cons_dir=="":
        print "Error: Need a bigwig file recording multiple alignment score."
    if mir_dir=="":
        print "Error: Need a gff3 file annotating miRNA."
