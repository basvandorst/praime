# Decisive comparator: re-map real zeros by p's explicit-formula wave.
#  ADD   (p | D): good-prime surrogate x solves x + dN_p(x) = y (y = real unfolded zeros), phase uniform (16 phases)
#  REMT  (p not| D): remove the TRUE wave: y = x + dN_p(x), dN_p(t) = (1/pi) arg(1 - chi(p) p^-1/2 e^{i t ln p})
#  REMR  (p not| D): same with a random phase (adds a second, independent wave)
#  ADDS  (p | D): ADD with the wave's sign flipped (dN -> -dN): tests the skewness (odd) part
import os, numpy as np, math, json, sys
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
SETS={'B':(['r134_B.npz'],100,600),'Dn':(['r136_Dn.npz'],100,600),'Fn':(['r138_Fn.npz'],100,600),
      'C':(['r134_C.npz'],600,1000),'Bp':(['r138_Bp.npz'],600,1000),
      'R120':(['r117_part0.npz','r117_part1.npz','r120_part0.npz','r120_part1.npz'],100,600)}
PR=(2,3,5,7,11,13); NPH=16
def kron(D,p):
    if D%p==0: return 0
    if p==2: return 1 if D%8 in (1,7) else -1
    return 1 if pow(D%p,(p-1)//2,p)==1 else -1
def unfold(t,q): return t/(2*math.pi)*np.log(q*t/(2*math.pi*math.e))
def tinv(x,q,t0):
    t=t0.copy()
    for _ in range(6): t=t+(x-unfold(t,q))*2*math.pi/np.log(q*t/(2*math.pi))
    return t
def wave(t,p,c,ph): return np.angle(1-c*p**-0.5*np.exp(1j*(t*math.log(p)+ph)))/math.pi
def st(x):
    s=np.diff(x); r=np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1]); return np.mean((s-s.mean())**2), r.mean()
def add_wave(t,q,p,ph,sgn=1):
    y=unfold(t,q); x=y.copy(); tt=t.copy()
    for _ in range(40):
        tt=tinv(x,q,tt); xn=y-sgn*wave(tt,p,1,ph)
        if np.max(abs(xn-x))<1e-12: x=xn; break
        x=xn
    return x
rng=np.random.default_rng(147)
out={}
names=sys.argv[1:] or list(SETS)
for name in names:
    files,lo,hi=SETS[name]; Z={}
    for f in files:
        z=np.load(LZ+f); Z.update({int(k[1:]):np.array(z[k]) for k in z.files})
    res={}
    for D,g in Z.items():
        q=abs(D); t=g[(g>=lo)&(g<=hi)]; x0=unfold(t,q); v0,r0=st(x0)
        rec={'v0':v0,'r0':r0}
        for p in PR:
            c=kron(D,p)
            if c==0:
                vv=[];rr=[];vs=[];rs=[]
                off=rng.uniform(0,2*math.pi)
                for j in range(NPH):
                    ph=off+2*math.pi*j/NPH
                    a,b=st(add_wave(t,q,p,ph)); vv.append(a); rr.append(b)
                    a,b=st(add_wave(t,q,p,ph,-1)); vs.append(a); rs.append(b)
                rec['add%d'%p]=(np.mean(vv)-v0,np.mean(rr)-r0); rec['adds%d'%p]=(np.mean(vs)-v0,np.mean(rs)-r0)
            else:
                xt=x0+wave(t,p,c,0.0); a,b=st(xt); rec['remt%d'%p]=(a-v0,b-r0)
                vv=[];rr=[]; off=rng.uniform(0,2*math.pi)
                for j in range(NPH):
                    xr=x0+wave(t,p,1,off+2*math.pi*j/NPH); a,b=st(xr); vv.append(a); rr.append(b)
                rec['remr%d'%p]=(np.mean(vv)-v0,np.mean(rr)-r0)
        res[str(D)]=rec
    out[name]=res
    print(name,'done',len(res),flush=True)
json.dump(out,open('remap_%s.json'%('_'.join(names)),'w'))
