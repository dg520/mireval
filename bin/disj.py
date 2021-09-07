def disjoint(l,res=[]):
    sl = []
    n = 1
    while n < len(l):
        p,q = l[n-1:n+1]
        if p[2] >= q[1] and p[0]==q[0]:
           l.remove(q)
           sl.append(q)
        else:
            n+=1
    res.append(l)
    if sl: 
        disjoint(sl,res)
    return res



l = [('chr1',100,300), ('chr1', 150, 350), ('chr1', 200, 400), ('chr1', 250, 450),('chr1', 500, 700), ('chr1',600, 800), ('chr1', 900, 1000)]
print disjoint(l,res=[])
l=[('chr1',100,300), ('chr1', 150, 350), ('chr1', 200, 400), ('chr1', 250, 450), ('chr1', 500, 700), ('chr1',600, 800), ('chr1', 900, 1000),('chr1',950,1100)]
print disjoint(l,res=[])
l=[('chr1',100,300), ('chr1', 150, 350), ('chr1', 200, 400), ('chr1', 250, 450), ('chr1', 500, 700), ('chr1',600, 800), ('chr1',650,880),('chr1', 900, 1000),('chr1',950,1100),('chr1',970,1200)]
print disjoint(l,res=[])
