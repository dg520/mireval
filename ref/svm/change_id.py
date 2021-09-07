from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sp_container=["homo_sapiens","mus_musculus"]

for sp in sp_container:
    print sp
    new=[]
    handle = open("/home/dadi/projects/web/genomes/"+sp+"/blastdb.fa", "rU")
    for record in SeqIO.parse(handle, "fasta"):
        t=record.description
        if ("mitochondrion" in t)==False:
            name="chr"+t.split(",")[0].split("chromosome ")[1]
        else:
            name="chrM"
        seq=record.seq
        content=SeqRecord(seq,name,"","")
        new.append(content)
        print name
    handle.close()

    output_handle = open("/home/dadi/projects/web/genomes/"+sp+"/blastdb_new.fa", "w")
    SeqIO.write(new, output_handle, "fasta")
    output_handle.close()
	