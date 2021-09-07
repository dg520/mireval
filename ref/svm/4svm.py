from __future__ import division
import linecache
import shlex, subprocess
import commands
import re



asb=["(((",".((","(.(","((.","..(",".(.","(..","..."]
nc=['a','u','c','g']
fea=[]
for a in nc:
    for b in asb:
        fea.append(a+b)

pf = raw_input("File: ") 
temp = open(pf,'r')
data=temp.readlines()
temp.close()


nop=raw_input("Negtive or Positive? ")
if nop=="+":
    nopw="+1"
if nop=="-":
    nopw="-1"

o=raw_input ("Output: ")

cmd="rm "+o
sta,out=commands.getstatusoutput(cmd)

temp = open(o,'w')
temp.close()



n=0
while n<len(data):
    seq =data[n+1][:-1]
    fold=data[n+2].split(" ")[0]
    fol=data[n+2].split(" ")[0]
    mfe=data[n+2].split(" ")[-1].split(")")[0].split("(")[-1]

    seq=seq.lower()
    
    i = 0 
    while i<len(fold):
	if fold[i]==")":
            fold=fold[:i]+"("+fold[(i+1):]
	i += 1
    

    mtx=[]
    m1x=[]
    i = 0
    while i<(len(seq)-2):
	ele=fold[i]+fold[i+1]+fold[i+2]
	mtx.append(ele)
	m1x.append(seq[i+1])
	i+=1

    train = open(o,'a')
    train.write(nopw)

    fx = 0
    for x in fea:
	i = 0
	val=0
	while i<len(mtx):
            if x[1:] == mtx[i] and x[0]==m1x[i]:
        	val += 1
            i += 1
	fx += 1
	fval = val/(len(seq)-2)
	train.write(" "+str(fx)+":"+str(fval))

    # write MFE
    fx=fx+1
    train.write(" "+str(fx)+":"+mfe)
    
    
    # write inter loop length
    
    fx=fx+1
    re1=re.compile("\(\.*\)")
    loop=re1.findall(fol)
    train.write(" "+str(fx)+":"+str((len(max(loop)))-2))
    
    
    #write bulge number
    fx=fx+1
   
    foldtemp=fol
    for lo in loop:
        foldtemp="".join(foldtemp.split(lo))          #eliminate inter loops from structure
    
    re3=re.compile("\([\(\.\)]*\)")
    h=re3.findall(foldtemp)
    newfold=h[0]                                      #eliminate bulges at two ends
    
    re2=re.compile("\.+")
    bulge=re2.findall(newfold)
    bulgen=len(bulge)
    train.write(" "+str(fx)+":"+str(bulgen))
    
    #write maximal bulge size
    fx=fx+1
    if bulgen != 0:
        train.write(" "+str(fx)+":"+str(len(max(bulge)))+"\n")
    else:
        train.write(" "+str(fx)+":"+"0\n")
    
    
    train.close()
    n+=3
       
