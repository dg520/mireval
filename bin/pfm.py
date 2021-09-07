f=open("matrix_only.txt","r")
d=f.readlines()
f.close()

i=0
while i<len(d):
    name=d[i][:-1].split(" ")[1]
    name=name.replace("/","or")
    f=open(name+".pfm","w")
    f.write("\t".join(d[i+1].split("[")[1].split("]")[0].split())+"\n")
    f.write("\t".join(d[i+2].split("[")[1].split("]")[0].split())+"\n")
    f.write("\t".join(d[i+3].split("[")[1].split("]")[0].split())+"\n")
    f.write("\t".join(d[i+4].split("[")[1].split("]")[0].split())+"\n")
    f.close()
    i+=5

