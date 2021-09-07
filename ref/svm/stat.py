f=open("res.txt","r")
d=f.readlines()
f.close()

dp=d[:715]
dn=d[715:]
len(dp)==len(dn)

tp=0
tn=0
fp=0
fn=0
for x in dp:
    t=float(x[:-1])
    if t>0:
        tp+=1
    else:
        fn+=1

for y in dn:
    t=float(y[:-1])
    if t<0:
        tn+=1
    else:
        fp+=1

print "True P", tp
print "False N", fn
print "True N", tn
print "False P", fp