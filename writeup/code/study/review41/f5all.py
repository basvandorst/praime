import cypari2, math, json
pari=cypari2.Pari(); pari.allocatemem(800*10**6); pari.default("realprecision",28)
pols=open('/home/bas/praime/writeup/code/study/lowzero/f5list.txt').read().split('\n')
rows=[]
for f in pols:
    P=pari(f); gal=pari.polgalois(P); assert int(gal[0])==20
    L=pari.lfundiv(pari.lfuncreate(P),pari.lfuncreate(1))
    N,_,mu=pari.lfunparams(L); N=int(N); mu=[int(m) for m in mu]
    assert N==abs(int(pari.nfdisc(P)))
    z=[float(x) for x in pari.lfunzeros(L,6)]
    v=float(pari.lfun(L,1)); d=float(pari.lfun(L,1,1))
    S=2*d/v+math.log(N/math.pi**4)+sum(float(pari.psi((1+m)/2)) for m in mu)
    rows.append(dict(f=f,N=N,mu=mu,g1=z[0],g2=z[1] if len(z)>1 else None,L1=v,S=S,irr=int(pari.polisirreducible(P))))
json.dump(rows,open('f5mine.json','w'))
print(len(rows), len(set(r['N'] for r in rows)), min(r['N'] for r in rows), max(r['N'] for r in rows))
