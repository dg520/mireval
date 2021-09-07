f=open("pol2.txt","r")
d=f.readlines()
f.close()


da=[]
dc=[]
dg=[]
dt=[]
line=[]
filen="name"
for x in d:
    t=x.split("\t")
    if t[0]!= filen:
        if dt!=[]:
            da=head+da
            line.append("\t".join(da)+"\n")
            line.append("\t".join(dc)+"\n")
            line.append("\t".join(dg)+"\n")
            line.append("\t".join(dt)+"\n")
            f.writelines(line)
            f.close()
            da=[]
            dc=[]
            dg=[]
            dt=[]
            line=[]
        f=open("/home/dadi/projects/web/jaspar/"+t[0]+".pfm","w")
        filen=t[0]
        head=[t[3].split(".")[0]]
    else:
        if t[1]=="A":
            da.append(t[3].split(".")[0])
        if t[1]=="C":           
            dc.append(t[3].split(".")[0])
        if t[1]=="G":
            dg.append(t[3].split(".")[0])
        if t[1]=="T":
            dt.append(t[3].split(".")[0])

da=head+da
line.append("\t".join(da)+"\n")
line.append("\t".join(dc)+"\n")
line.append("\t".join(dg)+"\n")
line.append("\t".join(dt)+"\n")
f.writelines(line)
f.close()