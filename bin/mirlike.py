from __future__ import division
import linecache
import shlex, subprocess
import commands
import re

a=raw_input("File: ")
temp = open(a,'r')
data=temp.readlines()
temp.close()
print len(data)

d=[]
n=0
while n<(len(data)-2):
    seq =data[n+1][:-1]
    fold=data[n+2].split(" ")[0]
    fol=data[n+2].split(" ")[0]
    mfe=data[n+2].split(" ")[-1].split(")")[0].split("(")[-1]

    seq=seq.lower()
    
    num=0
    for a in fold:
        if a =="(":
	    num+=1
    
    # write inter loop length
    
    re1=re.compile("\(\.*\)")
    loop=re1.findall(fol)
    if len(loop)==1 and num>=18 and float(mfe)<-15:
        d.append(data[n])
        d.append(data[n+1])
        d.append(data[n+2])
        
    n+=3

print len(d)

b=raw_input("Output: ")
f=open(b,"w")
for x in d:
    f.write(x)
f.close()
