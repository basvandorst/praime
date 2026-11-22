import cypari2, json
pari = cypari2.Pari(); pari.allocatemem(400*10**6)
out=[]
for a in range(3,5001):
    for D in (-a,a):
        if not pari.isfundamental(D): continue
        L = pari.lfuncreate(D)
        z = [float(x) for x in pari.lfunzeros(L, 30)]
        out.append(dict(D=D, z=z, L1=float(pari.lfun(L,1)), dL1=float(pari.lfun(L,1,1))))
json.dump(out, open('data.json','w'))
print(len(out), sum(1 for o in out if o['D']<0), sum(1 for o in out if o['D']>0))
