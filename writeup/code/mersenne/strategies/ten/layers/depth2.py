import math, sys, json, time
import cypari2  # needs cypari2 (PARI)
pari=cypari2.Pari(); pari.allocatemem(400*10**6)
W=int(sys.argv[1]); rows={}
t0=time.time()
for w in range(1,W+1):
    v=int(pari('polcyclo(%d,2)'%w)); 
    f=pari.factor(v); ps=[int(x) for x in f[0]]; es=[int(x) for x in f[1]]
    prim=[]
    for q,e in zip(ps,es):
        o=int(pari('znorder(Mod(2,%d))'%q)) if q>2 else 0
        if o==w: prim.append(q)
        assert all(pari.isprime(q) for q in ps)
    rows[w]=prim
print("primitive prime factors of 2^w - 1 (order of 2 exactly w), w = 1..%d, %.0f s"%(W,time.time()-t0))
tot=sum(len(v) for v in rows.values()); print("primes in the table:",tot,"; rows without one:",[w for w in rows if not rows[w]])
mers=[w for w in rows if len(rows[w])==1 and rows[w][0]==2**w-1]; print("rows whose only prime is 2^w - 1:",mers)
print("check of the floor: bits = log2(q+1) <= w for every q:",all(math.log2(q+1)<=w+1e-12 for w in rows for q in rows[w]),"; equality only for",[w for w in rows for q in rows[w] if q==2**w-1])
print("prime rows p: share of the bits of 2^p - 1 held by its largest prime factor")
for p in rows:
    if pari.isprime(p) and p>2:
        qs=sorted(rows[p]); print("  p=%3d: %d factor%s, bits %s, largest holds %.3f"%(p,len(qs),"" if len(qs)==1 else "s"," + ".join("%.1f"%math.log2(q) for q in qs),math.log2(qs[-1])/p))
json.dump({str(w):[str(q) for q in v] for w,v in rows.items()},open('depth_rows_'+str(W)+'.json','w'))
