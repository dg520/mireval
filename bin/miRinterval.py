f=open("mmu.gff3","r")
d=f.readlines()
f.close()

d1=[]
for x in d[11:]:
    if "miRNA_primary_transcript" in x:
        t=x.split("\t")
        chr="chr"+t[0]
        start=t[3]
        end=t[4]
        str=t[6]
        name=t[-1][:-1].split("=")[-1]
        d1.append(chr+"\t"+start+"\t"+end+"\t"+name+"\t0\t"+str+"\n")

f=open("mir_interval.txt","w")
f.writelines(d1)
f.close()
