import cypari2, math, random, sys, json, gmpy2
import numpy as np
pari=cypari2.Pari(); pari.allocatemem(800*10**6)
part=int(sys.argv[1]); npart=int(sys.argv[2]); N=int(sys.argv[3]) if len(sys.argv)>3 else 150
def g1(d):
    Lf=pari.lfuncreate(d); T=3*2*math.pi/math.log(abs(d)); z=[]
    while not z:
        z=[float(t) for t in pari.lfunzeros(Lf,T) if float(t)>1e-9]; T*=2
    return min(z)*math.log(abs(d))/(2*math.pi)
ds={}
for line in open('/home/bas/praime/writeup/code/study/polys/dist300.txt'):
    A,B,C,s,L=map(int,line.split())
    if L>=36:
        d=int(pari.coredisc(B*B-4*A*C)); ds.setdefault(d,[]).append(L)
P7=(2,3,5,7); P23=(2,3,5,7,11,13,17,19,23)
def kr(d,p): return int(gmpy2.kronecker(d,p))
def isfund(x):
    if x%4==1: return int(pari.issquarefree(x))==1
    if x%4==0:
        m=x//4
        return (m%4 in (2,3)) and int(pari.issquarefree(m))==1
    return False
out={}
rng=random.Random(9124+part)
for i,(d,runs) in enumerate(sorted(ds.items())):
    if i%npart!=part: continue
    gd=g1(d); sg=1 if d>0 else -1; lo,hi=int(abs(d)/1.5),int(abs(d)*1.5)
    pat7=[kr(d,p) for p in P7]; pat23=[kr(d,p) for p in P23]
    ref=[];m7=[];m23=[]
    tries=0
    cand=None
    if hi-lo<400000:
        cand=list(range(lo,hi+1)); rng.shuffle(cand)
    while (len(ref)<N or len(m7)<N or len(m23)<N):
        if cand is not None:
            if not cand: break
            x=sg*cand.pop()
        else:
            tries+=1
            if tries>1500000: break
            x=sg*rng.randint(lo,hi)
        if x==d or not isfund(x): continue
        k7=[kr(x,p) for p in P7]
        if len(ref)<N: ref.append(g1(x))
        if k7==pat7:
            if len(m7)<N: m7.append(g1(x))
            if len(m23)<N and [kr(x,p) for p in P23]==pat23: m23.append(g1(x))
    r=lambda a: float(np.mean(np.array(a)<gd)) if a else float('nan')
    out[d]=dict(runs=runs,g=gd,p_all=r(ref),p_m7=r(m7),p_m23=r(m23),med=[float(np.median(a)) if a else float('nan') for a in (ref,m7,m23)],n=[len(ref),len(m7),len(m23)])
    print(d,runs,'g %.4f'%gd,'pct all %.3f  match<=7 %.3f  match<=23 %.3f'%(r(ref),r(m7),r(m23)),'n',out[d]['n'],'medians',['%.3f'%x for x in out[d]['med']],flush=True)
json.dump({str(k):v for k,v in out.items()},open('r124mine_%d.json'%part,'w'))
