import random
from random import randint
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


a=["10","100","500","1000"]

handle=open("genomes/homo_sapiens/blastdb.fa","rU")
f=SeqIO.parse(handle,"fasta")
seq=[]
for x in f:
    seq.append(x.seq)
print "Done"

for y in a:
    n=1
    fa=[]
    for x in range(0,int(y)):
        id="id_"+str(n)
        no=random.sample(range(0,20),1)[0]
        start=randint(200000,len(seq[no])-200000)
        seq1=seq[no][start:start+91]
        record=SeqRecord(seq1,id,"","")
        fa.append(record)
        n+=1

    handle.close()

    output_handle = open("bin/ct"+y+".txt", "w")
    SeqIO.write(fa, output_handle, "fasta")
    output_handle.close()
    print y
