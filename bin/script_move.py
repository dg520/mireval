f=open("/srv/www/cgi-bin/mireval/mireval.py","r")
d=f.readlines()
f.close()

l=[]
for x in d[560:]:
    t="    "+x
    l.append(t)

d1=d[:560]+l



f=open("/srv/www/cgi-bin/mireval/mireval_new.py","w")
f.writelines(d1)
f.close()