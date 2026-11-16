import numpy as np, math, sys
from lib import sieve
N=20_000_000; s=sieve(N); P=s.astype(np.uint8); y0=int(round(N**0.4027)); sv=np.ones(N+1,dtype=bool); sv[:2]=False
for q in np.nonzero(s[:y0+1])[0]: sv[2*q::q]=False
BL=int(sys.argv[1]) if len(sys.argv)>1 else 12
nb=(N>>BL)+1; idx=np.arange(N+1)>>BL; bp=np.bincount(idx[y0+1:],weights=s[y0+1:],minlength=nb); bs=np.bincount(idx[y0+1:],weights=sv[y0+1:],minlength=nb); pb=np.where(bs>0,bp/np.maximum(bs,1),0)
rng=np.random.default_rng(11); ND=16; D=[]
for i in range(ND):
    d=(sv&(rng.random(N+1)<pb[idx])).astype(np.uint8); d[:y0+1]=P[:y0+1]; D.append(d)
n=np.arange(10**6,N-40000,dtype=np.int64); r=np.floor(np.sqrt(n)).astype(np.int64); r-=(r*r>n); r+=((r+1)**2<=n)
def test(h,label):
    m=n+h; np_=len(n)
    def R(L):
        u=L[n]; v=L[m]; pp=int((u&v).sum()); return pp,pp*np_/(int(u.sum())*int(v.sum()))
    pp,Rp=R(P); rd=[R(d) for d in D]; mean=np.mean([t[0] for t in rd]); Rd=np.mean([t[1] for t in rd]); sdR=np.std([t[1] for t in rd],ddof=1)
    z=(Rp-Rd)/(Rd*math.sqrt(1/pp+1/(ND*mean))); print("  %-22s pairs %7d  primes %.4f decoys %.4f (sd of one decoy %.4f)  z %+5.2f   z by decoy spread %+5.2f"%(label,pp,Rp,Rd,sdR,z,(Rp-Rd)/(sdR*math.sqrt(1+1/ND))))
print("blocks of %d numbers for the decoy density; n from 1e6 to 2e7"%(1<<BL))
for k in range(-4,5): test(2*r+2*k,"h = 2r %+d"%(2*k))
for k in range(-2,4): test(4*r+2*k,"h = 4r %+d"%(2*k))
for h in (6000,6002,6006,8190): test(np.int64(h),"h = %d (fixed)"%h)
rr=np.int64(1)*r
test(2*((r*7)//5),"h = 2*floor(1.4 r)")
test(2*r+2*(r%7),"h = 2r + 2(r mod 7)")
