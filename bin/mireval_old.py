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



# for i in *.gff3;do echo $i;head -n 10 $i | grep -m 1 -i "homo sapiens";done > toto
# less toto | grep -B 1 -i "Homo" 

# fetch seq from blastdb
# blastdbcmd -db blastdb -dbtype nucl -entry NC_000013 -strand minus -range 50623252-50623340 -outfmt %s

# fetch socre
# /home/shared/mireval/bigWigToBedGraph score.bw -chrom=chrX -start=60540 -end=60549 score.txt

# lift rfam and mir
# /home/shared/mireval/liftOver bta.gff3 -gff /home/shared/mireval/genomes/XXX/liftover/bosTau6ToBosTau7.over.chain new.gff3 unmapped.gff3


def conserve_score(root_name,result_no,chromosome,start_cor,end_cor,species,method):
    run=root_name+"/mireval/bin/bigWigToBedGraph "+root_name+"/mireval/genomes/"+species+method+"/score.bw -chrom=chr"+chromosome+" -start="+str(start_cor)+" -end="+str(end_cor)+" "+root_name+"/mireval/res/"+result_no+"/score_temp.txt"
    sta,out=commands.getstatusoutput(run)
    files=open(root_name+"/mireval/res/"+result_no+"/score_temp.txt","r")
    plot=files.readlines()
    files.close()
    if plot==[]:
        plot=[chromosome+"\t"+str(start_cor)+"\t"+str(end_cor)+"\t"+"0\n"]
    
    run="rm "+root_name+"/mireval/res/"+result_no+"/score_temp.txt"
    sta,out=commands.getstatusoutput(run)
    return plot

def score_zf(plots):
    zf=0
    for x in plots:
        tmpf=x[:-1].split("\t")
        zf=(int(tmpf[2])-int(tmpf[1])+1)*int(tmpf[3])+zf
    return zf

def neighbour_mir(root_name,chromosome,chromosome0,start_cor,end_cor,strand,species):
    files=open(root_name+"/mireval/genomes/"+species+"/mir/mir.gff3","r")
    mir_list1=files.readlines()
    files.close()
    mir_list="\n".join(mir_list1).split("#\n\n")[-1].split("\n\n")

    upstream=[]
    for premir in mir_list:
        if ("miRNA_primary_transcript" in premir)==False:
            mirchr=premir.split("\t")[0]
            if mirchr[0:3]!="chr":
                mirchr="chr"+mirchr
            if mirchr[-2:]==".1":
                mirchr=mirchr[:-2]
            new_chromosome="NA"
            if "_" in chromosome:
                new_chromosome="chr"+chromosome.split("_")[1]
            mirstart=int(premir.split("\t")[3])
            mirend=int(premir.split("\t")[4])
            mirname=premir.split("\t")[-1].split("Name=")[-1].split(";")[0]
            miracc=premir.split("\t")[-1].split("accession_number=")[-1].split(";")[0]
            mirst=premir.split("\t")[6]
            if mirchr==chromosome or mirchr==new_chromosome:
                if max(start_cor-1000,mirstart)<min(end_cor+1000,mirend):
                    upstream.append(mirname+"\t"+miracc+";"+chromosome0+":"+str(max(start_cor-1000,mirstart))+"-"+str(min(end_cor+1000,mirend))+":"+mirst+"\n")
     
    return upstream


def neighbour_nc(root_name,chromosome,chromosome0,start_cor,end_cor,strand,species):
    files=open(root_name+"/mireval/genomes/"+species+"/rfam/rfam_chr.txt","r")
    chr_list=files.readlines()
    files.close()

    neigh_rf=[]
    for rf in chr_list:
        rf_chr=rf.split("\t")[0]
        new_chromosome="NA"
        if "_" in chromosome:
            new_chromosome="chr"+chromosome.split("_")[1]
        if chromosome==rf_chr or rf_chr==new_chromosome:                                                       # or ("complete" in rf_chr)
            files=open(rf.split("\t")[1][:-1],"r")
            nc_list=files.readlines()
            files.close()
            for nc in nc_list:
                if nc[:2]!="##":
                    if max(start_cor-1000,int(nc.split("\t")[3]))<min(end_cor+1000,int(nc.split("\t")[4])):
                        rf_acc=nc.split("\t")[8].split(";")[1].split("=")[1]
                        rf_name="nc-"+nc.split("\t")[8].split(";")[2].split("=")[1]
                        rf_st=nc.split("\t")[6]
                        neigh_rf.append(rf_name+"\t"+rf_acc+";"+chromosome0+":"+str(max(start_cor-1000,int(nc.split("\t")[3])))+"-"+str(min(end_cor+1000,int(nc.split("\t")[4])))+":"+rf_st+"\n")

    return neigh_rf
  
    

def neighbour_trans(root_name,chromosome,chromosome0,start_cor,end_cor,strand,species):
    files=open(root_name+"/mireval/genomes/"+species+"/trans/trans.txt","r")
    trans_list1=files.readlines()
    files.close()
    trans_list=trans_list1[1:]

    neigh_trans=[]
    for trans in trans_list:
        transchr=trans.split("\t")[1]
        if transchr[0:3]!="chr":
            transchr="chr"+transchr
        transst=trans.split("\t")[2]
        transstart=int(trans.split("\t")[3])
        transend=int(trans.split("\t")[4])
        #transname=trans.split("\t")[5][:-1]
        transacc=trans.split("\t")[0]
        if transchr==chromosome and transst==strand:
            if max(start_cor-1000,transstart)<min(end_cor+1000,transend):
                neigh_trans.append(transacc+";"+chromosome0+":"+str(max(start_cor-1000,transstart))+"-"+str(min(end_cor+1000,transend))+":"+transst+"\n")

    return neigh_trans



def disect(root_name,seq_name,seq_chr,chromosome0,seq_start,seq_end,seq_st,sequence,mir_length_container):        
    my_dna = Seq(sequence, generic_dna)
    sequence1=my_dna.transcribe()
    seq_frag=[]
    if mir_length_container+10<len(str(sequence1)):
        a=0
        while a<len(sequence1)-100:
            win=mir_length_container
            t=str(sequence1[a:a+win])
            name=">"+seq_name+";"+chromosome0+":"+str(seq_start+a)+"-"+str(seq_start+a+win-1)+":"+seq_st+"\n"
            seq_frag.append(name)
            seq_frag.append(t+"\n")
            a+=10
    else:
        name=">"+seq_name+";"+chromosome0+":"+str(seq_start)+"-"+str(seq_end)+":"+seq_st+"\n"
        seq_frag.append(name)
        seq_frag.append(str(sequence1)+"\n")

    return seq_frag 


def polii(root_name,seq_name,chromosome,chromosome0,start_cor,end_cor,st,sequence):
    seq_polii=[]
    seq_temp=Seq(sequence,generic_dna)
    fast=SeqRecord(seq_temp,seq_name,'','')
    seq_polii.append(fast)
    
    return seq_polii


def sort_list(list1):
    tmp_list=[]
    for el in list1:
        t=el.split(";")
        new_el=(t[0]+";"+t[1].split(":")[0]+":",int(t[1].split(":")[1].split("-")[0]),"-"+t[1].split(":")[1].split("-")[1]+":"+t[1].split(":")[2])
        tmp_list.append(new_el)
        
    list1_sorted=sorted(tmp_list,key=itemgetter(1))
    
    final_list=[]
    for el in list1_sorted:
        t=el[0]+str(el[1])+el[2]
        final_list.append(t)

    return final_list

def disjoint(l,res=[]):
    sl = []
    n = 1
    while n < len(l):
        p,q = l[n-1:n+1]
        if int(p.split(";")[1].split(":")[1].split("-")[1]) >= int(q.split(";")[1].split(":")[1].split("-")[0]):
           l.remove(q)
           sl.append(q)
        else:
            n+=1
    res.append(l)
    if sl: 
        disjoint(sl,res)
    return res


def cir_hist(root_name, job_id, inq_id):
    plot_content=""
    for x in ["histo1.txt","histo2.txt"]:
        histhead="<plot>\ntype=histogram\n"
        histfile="file = "+root_name+"/mireval/res/"+job_id+"/"+inq_id+"/"+x+"\n"
        if x=="histo1.txt":
            r1="r1=0.85r\n"
            r0="r0=0.55r\n"
        else:
            r1="r1=0.5r\n"
            r0="r0=0.2r\n"
        histtail1="max=1\nmin=0\nstroke_type=bin\nextend_bin=no\n"
        histaxis="<axes>\n<axis>\nspacing=0.1r\ncolor=lgrey\nthickness=2\n</axis>\n</axes>\n"
        rules1="<rules>\n<rule>\ncondition=var(value)>=0.5\ncolor=dorange\nfill_color=orange\n</rule>\n"
        rules2="<rule>\ncondition=var(value)<0.5\ncolor=dgreen\nfill_color=green\n</rule>\n</rules>\n"
        histtail2="</plot>\n"
        plot_content=plot_content+histhead+histfile+r1+r0+histtail1+histaxis+rules1+rules2+histtail2

    return plot_content

def cir_ref(root_name, job_id, inq_id):
    plot_content=""
    heathead="<plot>\ntype=heatmap\n"
    heatfile="file = "+root_name+"/mireval/res/"+job_id+"/"+inq_id+"/ref.txt\n"
    color="color = spectral-9-div-rev\n"
    stroke_thickness="stroke_thickness = 1\n"
    stroke_color="stroke_color     = black\n"
    r0="r0= 0.925r\n"
    r1="r1= 0.925r-25p\n"
    plot_content=plot_content+heathead+heatfile+color+stroke_thickness+stroke_color+r0+r1+"</plot>\n"

    return plot_content

def cir_heatmap(root_name, job_id, inq_id):
    cmd="ls "+root_name+"/mireval/res/"+""+job_id+"/"+inq_id+"/heatmap*.txt"
    sta,out=commands.getstatusoutput(cmd)
    para=0
    plot_content=""
    for x in out.split("\n"):
        heathead="<plot>\ntype=heatmap\n"
        #heatfile="file = "+root_name+"/mireval/res/"+job_id+"/"+inq_id+"/"+x+"\n"
        heatfile="file = "+x+"\n"
        color="color = spectral-9-div-rev\n"
        stroke_thickness="stroke_thickness = 1\n"
        stroke_color="stroke_color     = black\n"
        r0="r0= "+str(1.05+0.075*para)+"r\n"
        r1="r1= "+str(1.05+0.075*para)+"r+25p\n"

        plot_content=plot_content+heathead+heatfile+color+stroke_thickness+stroke_color+r0+r1+"</plot>\n"
        para+=1

    return plot_content

def cir_text(root_name, job_id, inq_id):
    cmd="ls "+root_name+"/mireval/res/"+""+job_id+"/"+inq_id+"/name*.txt"
    sta,out=commands.getstatusoutput(cmd)
    para=0
    plot_content=""
    for x in out.split("\n"):
        texthead="<plot>\ntype=text\ncolor=black\n"
        #textfile="file = "+root_name+"/mireval/res/"+job_id+"/"+inq_id+"/"+x+"\n"
        textfile="file = "+x+"\n"
        r0="r0= "+str(1.05+0.075*para)+"r+30p\n"
        r1="r1= "+str(1.15+0.075*para)+"r+30p\n"
        #textlink="show_links=yes\nlink_dims=4p,4p,8p,4p,4p\nlink_thickness=2p\nlink_color=red\n"
        textlabel="label_size=24p\nlabel_parallel=yes\nlabel_font=condensed\npadding=0p\nrpadding=0p\n"
        plot_content=plot_content+texthead+textfile+r0+r1+textlabel+"</plot>\n"
        para+=1

    return plot_content


def cir(root_name, job_id, inq_id, chr_which):
    cir_content=""
    cir_head1="<colors>\n<<include "+root_name+"/mireval/circos-0.62/etc/colors.conf>>\n<<include "+root_name+"/mireval/circos-0.62/etc/brewer.conf>>\n</colors>\n"
    cir_head2="<fonts>\n<<include "+root_name+"/mireval/circos-0.62/etc/fonts.conf>>\n</fonts>\n"
    cir_head3="<<include "+root_name+"/mireval/conf/ideogram.conf>>\n<<include "+root_name+"/mireval/conf/ticks.conf>>\n"
    cir_head4="<image>\n<<include "+root_name+"/mireval/circos-0.62/etc/image.conf>>\n</image>\n"
    cir_head=cir_head1+cir_head2+cir_head3+cir_head4
    karyotype   = "karyotype = "+root_name+"/mireval/res/"+job_id+"/"+inq_id+"/karyo.txt\n"

    chromosomes_units = "chromosomes_units = 1\n"
    chromosomes_display_default = "chromosomes_display_default =  no\n"
   
    # which chromosomes to display
    chr_show=""
    for chr_id in chr_which:
        chr_show=chr_show+chr_id+";"
    chromosomes       = "chromosomes = "+   chr_show[:-1]   +"\n"

    cir_content=cir_content+cir_head+karyotype+chromosomes_units+chromosomes_display_default+chromosomes+"<plots>\n"
    cir_content=cir_content+cir_hist(root_name, job_id, inq_id)
    cir_content=cir_content+cir_ref(root_name, job_id, inq_id)
    cir_content=cir_content+cir_heatmap(root_name, job_id, inq_id)
    cir_content=cir_content+cir_text(root_name, job_id, inq_id)
    cir_content=cir_content+"</plots>\n<<include "+root_name+"/mireval/circos-0.62/etc/housekeeping.conf>>\n"


    f=open(root_name+"/mireval/res/"+job_id+"/"+inq_id+"/circos.conf","w")
    f.writelines(cir_content)
    f.close()
    cmd=root_name+"/mireval/circos-0.62/bin/circos -conf "+root_name+"/mireval/res/"+job_id+"/"+inq_id+"/circos.conf -outputdir "+root_name+"/mireval/res/"+job_id+"/"+inq_id+"/"
    sta,out=commands.getstatusoutput(cmd)
    
    #print out
    return 1




#<Main>

cgitb.enable()
print("Content-Type: text/html\n")
form=cgi.FieldStorage()

root="/home/shared"
f=open(root+"/mireval/record_log.txt","r")
record_step=f.readlines()
f.close()

res_no=str(int(record_step[0][:-1])+1)
cmd="mkdir "+root+"/mireval/res/"+res_no
sta,out=commands.getstatusoutput(cmd)

seq_file=raw_input("File: ")
sp=raw_input("Species: ")


if sp=="homo_sapiens":                        # GRCh37/hg19
    method_avbl=(1,1,1,1)
    mir_length=83    
if sp=="mus_musculus":                        # GRCm38/mm10
    method_avbl=(1,1,1,1)
    mir_length=87
if sp=="bos_taurus":                          # Btau_4.6.1/bosTau7
    method_avbl=(0,0,1,1)
    mir_length=80
if sp=="caenorhabditis_briggsae":             # cb3
    method_avbl=(0,0,1,1)
    mir_length=99
if sp=="caenorhabditis_elegans":              # WS220/ce10
    method_avbl=(1,0,1,1)
    mir_length=89
if sp=="canis_familiaris":                    # canFam2; not up-to-date
    method_avbl=(1,0,1,0)
    mir_length=69
if sp=="drosophila_melanogaster":             # BDGP_5.0/dm3
    method_avbl=(1,0,1,1)
    mir_length=94
if sp=="danio_rerio":                         # Zv9/danRer7
    method_avbl=(1,1,1,0)
    mir_length=93
if sp=="gallus_gallus":                       # galGal4/ICGSC Gallus_gallus-4.0
    method_avbl=(0,0,1,1)
    mir_length=93
if sp=="xenopus_tropicalis":                  # xenTro3/JGI_4.2
    method_avbl=(1,0,1,0)
    mir_length=81
if sp=="rattus_norvegicus":                   # rn4/RGSC_3.4; not up-to-date
    method_avbl=(1,0,1,1)
    mir_length=91
if sp=="pan_troglodytes":                     # CSAC_2.1.4/panTro4
    method_avbl=(0,0,1,1)
    mir_length=89
if sp=="macaca_mulatta":                      # rheMac2/MMUL_1.0; not up-to-date
    method_avbl=(0,0,1,1)
    mir_length=89
if sp=="monodelphis_domestica":               # MonDom5
    method_avbl=(1,0,1,1)
    mir_length=80
if sp=="fugu_rubripes":                       # FUGU5/fr3
    method_avbl=(1,0,1,0)
    mir_length=81

if sp=="pongo_albelii":                       # ponAbe2/WUSTL_2.0.2
    method_avbl=(1,0,0,0)
    mir_length=89
if sp=="gorilla_gorilla":                     # gorGor3.1/gorGor3
    method_avbl=(1,0,0,0)
    mir_length=89
if sp=="branchiostoma_floridae":              # braFlo1/JGI_1.0
    method_avbl=(1,0,0,0)
    mir_length=93
if sp=="petromyzon_marinus":                  # WUGSC 7.0/petMar2
    method_avbl=(1,0,0,0)
    mir_length=93
if sp=="oryzias_latipes":                     # oryLat2/MEDAKA1
    method_avbl=(1,0,0,0)
    mir_length=93
if sp=="geospiza_fortis":                     # GeoFor_1.0/geoFor1
    method_avbl=(1,1,0,0)
    mir_length=93
if sp=="ornithorhynchus_anatinus":            # ornAna1/WUSTL_5.0.1
    method_avbl=(1,0,0,1)
    mir_length=87
if sp=="gasterosteus_aculeatus":              # gasAcu1/BI_1.0
    method_avbl=(1,0,0,0)
    mir_length=93
if sp=="saccharomyces_cerevisiae":            # sacCer3/Apr2011
    method_avbl=(1,0,0,1)
    mir_length=93

if sp=="sus_scrofa":                          # SGSC_Sscrofa10.2/susScr3
    method_avbl=(0,0,0,1)
    mir_length=80
if sp=="felis_catus":                         # Felis_catus_6.2/felCat5
    method_avbl=(0,0,0,0)
    mir_length=69
if sp=="oryctolagus_cuniculus":               # Broad/oryCun2
    method_avbl=(0,0,0,1)
    mir_length=89
if sp=="ovis_aries":                          # Ovis_aries_1.0/oviAri1
    method_avbl=(0,0,0,1)
    mir_length=80
if sp=="equus_caballus":                      # equCab2
    method_avbl=(0,0,0,1)
    mir_length=80
if sp=="macropus_eugenii":                    # TWGS_Meug_1.1/macEug2
    method_avbl=(0,0,0,0)
    mir_length=69
if sp=="sarcophilus_harrisii":                # Devil_ref v7.0/sarHar1
    method_avbl=(0,0,0,0)
    mir_length=89

# command to run blast
# makeblastdb -in blastdb.fa -dbtype nucl -parse_seqids -out blastdb

# fetch chromosome sizes
# /home/shared/mireval/fetchChromSizes mm10 > chrom.sizes

# wig to bigWIG; not in program
# zcat *.gz | gzip -c > score.wigFix.gz
# cd..
# /home/shared/mireval/wigToBigWig cons/score.wigFix.gz chrom.sizes cons/score.bw


hit_thres=0.95


cmd="blastn -db "+root+"/mireval/genomes/"+sp+"/blastdb -query "+root+"/mireval/"+seq_file+" -outfmt 5 -out "+root+"/mireval/res/"+res_no+"/"+"temp.txt"
sta,out=commands.getstatusoutput(cmd)
if out!="":
    print("<html><body>Error: Your input data is not in FASTA format. Please use the correct Format and try again.")

else:
 
blast_records = NCBIXML.parse(open(root+"/mireval/res/"+res_no+"/temp.txt"))


dblast={}
dfrag=[]
dpolii=[]
id_no=0
for blast_record in blast_records:
    blast_no=0
    name=str(blast_record.query)

    dblast[name]=[]

    hit_len=blast_record.query_letters


    chrom_file=open(root+"/mireval/genomes/"+sp+"/chrom.sizes")
    chrom_size=chrom_file.readlines()
    chrom_file.close()

    hn=0
    for alignment in blast_record.alignments:
        dict={"inq":"", "intra_id":"", "no":0, "chr":"", "chr0":"", "start":0, "end":0, "start_outer":0, "end_outer":0, "strand":"" ,"alig":[],"plot1":[], "plot2":[], "promoter":[], "nearby":[], "frag":[], "pred":[], "trans":[],"nc":[],"disj":[]}
        hsp=alignment.hsps[0]
        if hsp.identities<hit_len*hit_thres:
            if hn==0:
                dblast[name].append(dict)
            break
        if hsp.identities>=hit_len*hit_thres:
            hn=1
            blast_no+=1
            id_no+=1
            dict["intra_id"]="id_"+str(id_no)
            dict["no"]=blast_no
            blast_name=name+"_"+str(blast_no)
            query_start=int(hsp.query_start)
            query_end=int(hsp.query_end)
            entry=str(alignment.accession)
            start=min(int(hsp.sbjct_start),int(hsp.sbjct_end))
            end=max(int(hsp.sbjct_start),int(hsp.sbjct_end))
            dict["start"]=start
            dict["end"]=end
            if hsp.frame==(1,1):
                st="+"
                st_name="plus"
            else:
                st="-"
                st_name="minus"
            dict["strand"]=st
            chr=str((alignment.title).split("chromosome ")[-1].split(",")[0]).split(" No definition line")[0]
            dict["chr0"]=chr
            if chr[0:3]!="chr":
                chr1="chr"+chr
            else:
                chr1=chr
            dict["chr"]=chr1
            seq=str(hsp.sbjct).upper()
            
            for x in chrom_size:
                if chr==x.split("\t")[0]:
                    chrom_end=int(x.split("\t")[1][:-1])
                    break
            end_o=min(end+1000, chrom_end)
            start_o=max(start-1000,0)
            dict["start_outer"]=start_o
            dict["end_outer"]=end_o

            cmd="blastdbcmd -db "+root+"/mireval/genomes/"+sp+"/blastdb -dbtype nucl -entry "+entry+" -strand "+st_name+" -range "+str(start_o)+"-"+str(end_o)+" -outfmt %s"
            sta,out=commands.getstatusoutput(cmd)
            seq_around=out[:-1]
            dict["alig"].append("Query "+str(query_start)+"-"+str(query_end)+"\n"+chr+":"+str(start)+"-"+str(end)+":"+st+"\n"+"Strand=Plus/"+st_name+"\n"+str(hsp.query)+"\n"+str(hsp.match)+"\n"+str(hsp.sbjct)+"\n")
            if method_avbl[0]==1:          
                dict["plot1"]=conserve_score(root,res_no,chr,start,end,sp,"/cons")                                       # conservation res
            if method_avbl[1]==1:  
                dict["plot2"]=conserve_score(root,res_no,chr,start,end,sp,"/shadow")                                     # shadow res
            if method_avbl[2]==1:
                dict["nearby"]=neighbour_mir(root,chr1,chr,start,end,st,sp)
            if method_avbl[3]==1:
                dict["nc"]=neighbour_nc(root,chr1,chr,start,end,st,sp)
            dict["trans"]=neighbour_trans(root,chr1,chr,start,end,st,sp)
            
            dblast[name].append(dict)
            dfrag=dfrag+disect(root,"id_"+str(id_no),chr1,chr,start,end,st,seq,mir_length)
            dpolii=dpolii+polii(root,"id_"+str(id_no),chr1,chr,start,end,st,seq_around)
    


# Find + prediction by HeteroMirPred

dfrag1=[]
d=[]
m=0
for x in dfrag:
    if m%2==0:
        d.append(x[:-1])
        dfrag1.append(">name"+str(m)+"\n")
    else:
        dfrag1.append(x)
    m+=1


f=open(root+"/mireval/bin/res_frag_"+res_no+".fasta","w")
f.writelines(dfrag1)
f.close()

cmd="perl HeteroMirPred.pl res_frag_"+res_no+".fasta"
sta,out=commands.getstatusoutput(cmd)


preds=out.split("##\n\n")[-1].split("\n")


positive=[]
pred_no=1
fasta_no=0
for pred in preds:
    if pred.split(":")[-1].split(" ")[0]=="miRNA":
        positive.append((fasta_no,pred_no))
        pred_no+=1
    fasta_no+=1


struct2nd_rec=out.split("zdata\n")[1:]



for x in positive:
    for y in dblast:
        for z in dblast[y]:
            if z["intra_id"]==d[x[0]].split(";")[0][1:]:
                pred_name="pred"+str(x[1])+";"+d[x[0]].split(";")[1]
                z["pred"].append(pred_name+"\n")
                z["frag"].append(pred_name+"\n")
                z["frag"].append(struct2nd_rec[x[0]].split("\n")[0]+"\n")
                z["frag"].append(struct2nd_rec[x[0]].split("\n")[1]+"\n")
                break


cmd="rm "+root+"/mireval/bin/res_frag_"+res_no+".fasta*"
sta,out=commands.getstatusoutput(cmd)


# Promoter Finder
f=open(root+"/mireval/res/"+res_no+"/res_4fimo_"+res_no+".txt","w")
SeqIO.write(dpolii, f, "fasta")
f.close()

cmd=root+"/mireval/meme/bin/fimo --output-pthresh 1e-5 --oc "+root+"/mireval/res/"+res_no+"/fimo_out/ "+root+"/mireval/jaspar.meme "+root+"/mireval/res/"+res_no+"/res_4fimo_"+res_no+".txt" 
sta,out=commands.getstatusoutput(cmd)


f=open(root+"/mireval/res/"+res_no+"/fimo_out/fimo.txt","r")
d=f.readlines()
f.close()


for x in d[1:]:
    for y in dblast:
        for z in dblast[y]:
            if z["intra_id"] == x.split("\t")[1]:
                pro_name=x.split("\t")[0]
                pro_num1=z["start_outer"]+int(x.split("\t")[2])-1
                pro_num2=z["start_outer"]+int(x.split("\t")[3])-1
                pro_start=min(pro_num1,pro_num2)
                pro_end=max(pro_num1,pro_num2)
                z["promoter"].append("promoter-"+pro_name+";"+z["chr"]+":"+str(pro_start)+"-"+str(pro_end)+":"+z["strand"]+"\n")
                break



# create circos reference files, html file and text result file

print("<html><head></head><body><table border=1 align=\"center\"><tr valign=\"middle\"><h2 align=\"center\">Inquiry ID: 10001</h2>")
print("Download<br><br><br>")
print("<b>Pred:</b> In red if miREval finds a predicted miRNA<br>")
print("<b>miRB:</b> In red if miREval finds a miRBase miRNA<br>")
print("<b>rFam:</b> In red if miREval finds a fFam non-coding RNA<br>")
print("<b>Cons:</b> In red if the genomic conservation score is avaibale<br>")
print("<b>Shad:</b> In red if the phylogenetic shadowing is avaibale<br><br>")
print("<table border=\"1\">")
print("<tr><th>Name</th><th>miRB</th><th>rFam</th><th>Pred</th><th>Cons</th><th>Shad</th></tr>")

text_res=open(root+"/mireval/res/"+res_no+"/res_txt.txt","w")
for x in dblast:
    kary=[]
    ref=[]
    transcripts=[]
    heatmap=[]
    name=[]
    hist1=[]
    hist2=[]
    leng=[]
    all_neigh=[]
    amount_neigh=[]
    all_pred=[]
    all_nc=[]
    amount_nc=[]
    chr_display=[]
    all_alig=[]
    all_promoter=[]
    disjo=[]


    cmd="mkdir "+"\""+root+"/mireval/res/"+res_no+"/"+x+"\""
    sta,out=commands.getstatusoutput(cmd)

    print("<tr><td height=20><a href=\"res/1.txt\">"+x+"</a></td>")
    ini=0
    for y in dblast[x]:
        if y["chr"]=="":
            htm.write("<html><head><title>Result of "+x+"</title></head>"+"<body><h2 align=\"center\">Result of "+x+"</h2>"+"<b>Inquiry ID:</b><br>"+res_no+"<br><br><br>"+"<b>Inquiry:</b><br>")
            htm.write("<b>Not Mappable</b><br>")
            htm.close()
            text_res.write(x+" start\n")
            text_res.write("Not mappable\n")
            text_res.write(x+" end\n\n\n")
            break
        else:
            ini=1
            if (y["chr"] in chr_display)==False:
                chr_display.append(y["chr"])
            kary.append("chr"+"\t"+"-"+"\t"+y["chr"]+"\t"+y["chr"]+"\t"+str(y["start_outer"])+"\t"+str(y["end_outer"])+"\t"+"grey"+"\t"+y["strand"]+"\n")
            ref.append(y["chr"]+"\t"+str(y["start"])+"\t"+str(y["end"])+"\t100\t"+"color=yellow"+"\t"+y["strand"]+"\n")
            transcripts=transcripts+y["trans"]
            hist1=hist1+y["plot1"]
            hist2=hist2+y["plot2"] 
            a=y["nearby"]+y["promoter"]+y["pred"]+y["nc"]+y["trans"]
            aa=sort_list(a)
            disjo=disjoint(aa,res=[])
            heatmap.append(disjo)
            leng.append(len(disjo)) 
            all_neigh= all_neigh+y["nearby"]
            amount_neigh.append(len(y["nearby"]))
            all_pred=all_pred+y["frag"]
            all_nc=all_nc+y["nc"]
            amount_nc.append(len(y["nc"]))
            all_alig=all_alig+y["alig"]
            all_promoter=all_promoter+y["promoter"]


    if ini==1:
        text_res.write(x+" start\n")       
        htm=open(root+"/mireval/res/"+res_no+"/"+x+"/res.html","w")
        htm.write("<html><head><title>Result of "+x+"</title></head>"+"<body><h2 align=\"center\">Result of "+x+"</h2>"+"<b>Inquiry ID:</b><br>"+res_no+"<br><br><br>"+"<b>Inquiry:</b><br>")
        text_res.write("Inquiry ID: "+res_no+"\n")
        text_res.write("Inquiry:\n")
        for blast_info in all_alig:
            t0=blast_info.split("\n")
            htm.write(t0[0]+"<br>"+t0[1]+"<br>"+t0[2]+"<br>"+t0[3]+"<br>"+t0[4]+"<br>"+t0[5]+"<br><br>")
            text_res.write(blast_info)
        htm.write("<br><br><b>Transcripts:</b><br>")
        text_res.write("Transcripts:\n")
        
        if transcripts!=[]:
            for trans_info in transcripts:
                name_acc=trans_info.split(";")[0]
                t0=trans_info.split(";")[1][:-1]
                htm.write("<b>"+name_acc+"</b>:"+t0+"<br>")
                text_res.write(name_acc+"\t"+t0+"\n")
        else:
            htm.write("None<br>")
            text_res.write("None\n")

        htm.write("<br><br><b>Promoter:</b><br>")
        text_res.write("Promoter:\n")

        if all_promoter!=[]:
            for promoter_info in all_promoter:
                name=promoter_info.split(";")[0]
                t0=promoter_info.split(";")[1][:-1]
                htm.write("<b>"+name+"</b>:"+t0+"<br>")
                text_res.write(name+"\t"+t0+"\n")
        else:
            htm.write("None<br>")
            text_res.write("None\n")

        htm.write("<br><br><b>miRNA cluster from miRBase(ver_19):</b><br>")
        text_res.write("miRNA cluster from miRBase(ver_19):\n")
        
        if all_neigh!=[]:
            for mir_info in all_neigh:
                name_acc=mir_info.split(";")[0]
                name=name_acc.split("\t")[0]
                acc=name_acc.split("\t")[1]
                t0=mir_info.split(";")[1][:-1]
                htm.write("<a href=http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc="+acc+">"+acc+"</a>:<b>"+name+"</b>:"+t0+"<br>")
                text_res.write(acc+"\t"+name+"\t"+t0+"\n")
            if max(amount_neigh)>5:
                print("<td bgcolor=\"#B40404\"></td>")
            if max(amount_neigh)>3 and max(amount_neigh)<=5:
                print("<td bgcolor=\"#DF0101\"></td>")
            if max(amount_neigh)>1 and max(amount_neigh)<=3:
                print("<td bgcolor=\"#FA5858\"></td>")
            if max(amount_neigh)==1:
                print("<td bgcolor=\"#F5A9A9\"></td>")
        else:
            if method_avbl[2]==0:
                htm.write("NA<br>")
                text_res.write("NA\n")
                print("<td bgcolor=\"#F6CECE\"></td>")
            else:
                htm.write("None<br>")
                text_res.write("None\n")
                print("<td></td>")

        htm.write("<br><br><b>Non-coding RNA cluster from Rfam(ver_11.0):</b><br>")
        text_res.write("Non-coding RNA clustering from Rfam(ver_11.0):\n")

        if all_nc!=[]:
            for nc_info in all_nc:
                name_acc=nc_info.split(";")[0]
                name=name_acc.split("\t")[0]
                acc=name_acc.split("\t")[1]
                t0=nc_info.split(";")[1][:-1]
                htm.write("<a href=hhttp://rfam.sanger.ac.uk/family/"+acc+">"+acc+"</a>:<b>"+name+"</b>:"+t0+"<br>")
                text_res.write(acc+"\t"+name+"\t"+t0+"\n")
            if max(amount_nc)>5:
                print("<td bgcolor=\"#B40404\"></td>")
            if max(amount_nc)>3 and max(amount_nc)<=5:
                print("<td bgcolor=\"#DF0101\"></td>")
            if max(amount_nc)>1 and max(amount_nc)<=3:
                print("<td bgcolor=\"#FA5858\"></td>")
            if max(amount_nc)==1:
                print("<td bgcolor=\"#F5A9A9\"></td>")
        else:
            if method_avbl[3]==0:
                htm.write("NA<br>")
                text_res.write("NA\n")
                print("<td bgcolor=\"#F6CECE\"></td>")
            else:
                htm.write("None<br>")
                text_res.write("None\n")
                print("<td></td>")

        htm.write("<br><br><b>miRNA prediction:</b><br>")
        text_res.write("miRNA predictions:\n")
        
        # /home/shared/mireval/bigWigToBedGraph score.bw -chrom=chrX -start=60540 -end=60549 score.txt
        # conserve_score(root_name,result_no,chromosome,start_cor,end_cor,species,method)
        if all_pred!=[]:
            print("<td bgcolor=\"#B40404\"></td>")
            m=0
            cons_val=[]
            shad_val=[]
            for pred_info in all_pred:
                if m%3==0:
                    t0=pred_info.split(";")
                    htm.write(">"+t0[0]+":"+t0[1][:-1]+"<br>")
                    text_res.write(">"+t0[0]+":"+t0[1][:-1]+"\n")
                    if method_avbl[0]==1:
                        pred_q_chr=t0[1].split(":")[0]
                        pred_q_start=int(t0[1].split(":")[1].split("-")[0])
                        pred_q_end=int(t0[1].split(":")[1].split("-")[1])
                        pred_q_cons=conserve_score(root,res_no,pred_q_chr,pred_q_start,pred_q_end,sp,"/cons")
                        consv=score_zf(pred_q_cons)/(pred_q_end-pred_q_start+1)
                    else:
                        pred_q_cons="NA"
                    cons_val.append(consv)
                    if method_avbl[1]==1:
                        pred_q_chr=t0[1].split(":")[0]
                        pred_q_start=t0[1].split(":")[1].split("-")[0]
                        pred_q_end=t0[1].split(":")[1].split("-")[1]
                        pred_q_shadow=conserve_score(root,res_no,pred_q_chr,pred_q_start,pred_q_end,sp,"/shadow")
                        shadv=score_zf(pred_q_shadow)/(pred_q_end-pred_q_start+1)
                    else:
                        pred_q_shadow="NA"
                    shad_val.append(shadv)
                if m%3==1:
                    pred_seq=pred_info[:-1]
                    htm.write(pred_info[:-1]+"<br>")
                    text_res.write(pred_seq+"\n")
                if m%3==2:
                    pred_struct=pred_info[:-1].split(" ")[0]
                    htm.write(pred_info[:-1]+"<br>")
                    text_res.write(pred_info[:-1]+"\n")
                    htm.write("<applet  code=\"VARNA.class\" codebase=\""+root+"/mireval/varna\" "+ "archive=\"VARNAv3-9.jar\" width=\"720\" height=\"600\">")
                    htm.write("<param name=\"title\" value=\""+t0[0]+":"+t0[1][:-1]+"\"/><param name=\"titleSize\" value=\"18\" />")
                    htm.write("<param name=\"sequenceDBN\" value=\""+pred_seq+"\"/>")
                    htm.write("<param name=\"structureDBN\" value=\""+pred_struct+"\"/></applet><br><br>")
                m+=1
            if max(cons_val)>0.9:
                print("<td bgcolor=\"#B40404\"></td>")
            if max(cons_val)>0.7 and max(cons_val)<=0.9:
                print("<td bgcolor=\"#DF0101\"></td>")
            if max(cons_val)>0.5 and max(cons_val)<=0.7:
                print("<td bgcolor=\"#FA5858\"></td>")
            if max(cons_val)<=0.5:
                print("<td bgcolor=\"#F5A9A9\"></td>")
            if max(cons_val)=="NA":
                print("<td></td>")
            if max(shad_val)>0.9:
                print("<td bgcolor=\"#B40404\"></td>")
            if max(shad_val)>0.7 and max(shad_val)<=0.9:
                print("<td bgcolor=\"#DF0101\"></td>")
            if max(shad_val)>0.5 and max(shad_val)<=0.7:
                print("<td bgcolor=\"#FA5858\"></td>")
            if max(shad_val)<=0.5:
                print("<td bgcolor=\"#F5A9A9\"></td>")
            if max(shad_val)=="NA":
                print("<td></td>")
        else:
            htm.write("None<br>")
            text_res.write("None\n")
            print("<td></td>")
            print("<td></td>")
            print("<td></td>")
        htm.write("<br><br><b>Phylogenic Conservation Score:</b><br>")
        text_res.write("Phylogenic Conservation Score:\n")


        if hist1!=[]:
            for con_info in hist1:
                t0=con_info.split("\t")
                htm.write(t0[0]+":"+t0[1]+"-"+t0[2]+":"+t0[3][:-1]+"<br>")
                text_res.write(t0[0]+":"+t0[1]+"-"+t0[2]+":"+t0[3][:-1]+"\n")
        else:
            htm.write("NA<br>")
            text_res.write("NA\n")
        htm.write("<br><br><b>Phylogenic Shadowing Score:</b><br>")
        text_res.write("Phylogenic Shadowing Score:\n")

        if hist2!=[]:
            for shadow_info in hist2:
                t0=shadow_info.split("\t")
                htm.write(t0[0]+":"+t0[1]+"-"+t0[2]+":"+t0[3][:-1]+"<br>")
                text_res.write(t0[0]+":"+t0[1]+"-"+t0[2]+":"+t0[3][:-1]+"\n")
        else:
            htm.write("NA<br>")
            text_res.write("NA\n")

        htm.write("<br><br><div align=\"center\">"+"<img align=\"center\" src=\"circos.png\" width=\"2000\" height=\"2000\">"+"</div></body></html>")

#
        f=open(root+"/mireval/res/"+res_no+"/"+x+"/karyo.txt","w")
        f.writelines(kary)
        f.close()
        f=open(root+"/mireval/res/"+res_no+"/"+x+"/ref.txt","w")
        f.writelines(ref)
        f.close()
        f=open(root+"/mireval/res/"+res_no+"/"+x+"/histo1.txt","w")
        f.writelines(hist1)
        f.close()
        f=open(root+"/mireval/res/"+res_no+"/"+x+"/histo2.txt","w")
        f.writelines(hist2)
        f.close()       
        for i in range(0,max(leng)):
            if i>0:
                f.close()
                fi.close()
            f=open(root+"/mireval/res/"+res_no+"/"+x+"/heatmap"+str(i+1)+".txt","w")
            fi=open(root+"/mireval/res/"+res_no+"/"+x+"/name"+str(i+1)+".txt","w")
            for each in heatmap:
                if i< len(each):
                    for line in each[i]:
                        type_info=line.split(";")[0].split("\t")[0]
                        if ("promoter-" in type_info) or ("pred" in type_info) or ("miR-" in type_info) or ("nc-" in type_info):
                            if "promoter-" in type_info:
                                type_info= type_info.split("-")[1]
                                color_info="color=green"
                            if "pred" in type_info:
                                color_info="color=vdred"
                            if "miR-" in type_info:
                                type_info=type_info[4:]
                                color_info="color=blue"
                            if "nc-" in type_info:
                                color_info="color=orange"
                        else:
                            color_info="color=purple"
                        t1=line.split(";")[1].split(":")
                        chr_info=t1[0]
                        start_info=t1[1].split("-")[0]
                        end_info=t1[1].split("-")[1]
                        st_info=t1[2][0]
                        f.write(chr_info+"\t"+start_info+"\t"+end_info+"\t"+"100\t"+color_info+"\t"+st_info+"\n")
                        fi.write(chr_info+"\t"+start_info+"\t"+end_info+"\t"+type_info+"\n")
            f.close()
            fi.close()
        circos_create=cir(root, res_no, x,chr_display)
        text_res.write(x+" end\n\n\n")

text_res.close()

print("</table>")



# record current user-ID
f=open(root+"/mireval/record_log.txt","w")
f.writelines([res_no+"\n"])
f.close()

print ("</body></html>")


# circos-0.62/bin/circos -conf res/10001/circos.conf -outputdir res/

# function to show histgram                            (OK)
# on circos, certain chromosomes to show               (OK)
# rotate the name of each circos bar                   (OK)
# on web, java of 2nd structure                        (OK)
# text result for download                             (OK)
# core pipeline debugging                              (OK)
# all promoter motifs establishment                    (OK)
# other species                                        (OK)
# queue system                                         (OK)
# temprary page to show waiting and dealing            ()
# visit count                                          ()
