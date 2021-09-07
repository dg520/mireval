def fastx(a):
    new=[]
    seq=""
    i=0
    while i<len(a):
        if a[i][0]==">":
            new.append(seq+"\n")
            new.append(a[i])
            seq=""
            fg=1
        else:
            seq=seq+a[i][:-2]
        i+=1
    new.append(seq+"\n")
    return new[1:]


f=open("rn45s.txt","r")
data=f.readlines()
f.close()
d=fastx(data)
f=open("rn45s.fastq","w")
f.writelines(d)
f.close()


