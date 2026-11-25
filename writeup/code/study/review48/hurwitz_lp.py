# independent L'(rho, chi_D) by Hurwitz zeta (PARI zetahurwitz), not lfun
import cypari2, json, sys, time, os
pari=cypari2.Pari(); pari.set_real_precision(30); pari.allocatemem(200*10**6)
jobs=json.load(open(sys.argv[1])); outf=sys.argv[2]; tmax=float(sys.argv[3])
done=set()
if os.path.exists(outf):
    for l in open(outf): done.add(tuple(json.loads(l)[:3]))
pari('Lv(D,s)=my(q=abs(D),A=0,B=0,c);for(a=1,q,c=kronecker(D,a);if(c,A+=c*zetahurwitz(s,a/q);B+=c*zetahurwitz(s,a/q,1)));A*=q^(-s);B*=q^(-s);[A,B-log(q)*A]')
t0=time.time()
with open(outf,'a') as fo:
    for s,D,g,y in jobs:
        if (s,D,g) in done: continue
        if time.time()-t0>tmax: break
        r=pari('Lv(%d,1/2+I*%.15f)'%(D,g))
        L=complex(r[0]); Lp=complex(r[1])
        fo.write(json.dumps([s,D,g,y,abs(L),abs(Lp)])+'\n'); fo.flush()
print('done',len(done),'->',sum(1 for _ in open(outf)),'of',len(jobs), 'time',time.time()-t0)
