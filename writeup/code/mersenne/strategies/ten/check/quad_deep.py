import numpy as np, math, sys
from lib import sieve
N=20_000_000; s=sieve(N); P=s.astype(np.uint8); y0=int(round(N**0.4027)); sv=np.ones(N+1,dtype=bool); sv[:2]=False
for q in np.nonzero(s[:y0+1])[0]: sv[2*q::q]=False
nb=(N>>12)+1; idx=np.arange(N+1)>>12; bp=np.bincount(idx[y0+1:],weights=s[y0+1:],minlength=nb); bs=np.bincount(idx[y0+1:],weights=sv[y0+1:],minlength=nb); pb=np.where(bs>0,bp/np.maximum(bs,1),0)
rng=np.random.default_rng(7); D=[]
for i in range(8):
    d=(sv&(rng.random(N+1)<pb[idx])).astype(np.uint8); d[:y0+1]=P[:y0+1]; D.append(d)
print("N=%d, sieve to %d, primes %d, decoy sizes %s"%(N,y0,int(P.sum()),[int(d.sum()) for d in D]))
n=np.arange(1,N+1,dtype=np.int64); r=np.floor(np.sqrt(n)).astype(np.int64); r-=(r*r>n); r+=((r+1)**2<=n); x=n-r*r
def run(lo,hi,label):
    sel=(n>lo)&(n<=hi); nn=n[sel]; rr=r[sel]; xx=x[sel]; zz=[]; out=[]
    for dy in range(0,4):
        for dx in range(-3,4):
            if dy==0 and dx<=0: continue
            r2=rr+dy; x2=xx+dx; ok=(x2>=0)&(x2<=2*r2); m=r2*r2+x2; ok&=(m<=hi)&(m>lo); a_=nn[ok]; m_=m[ok]; np_=len(a_)
            def R(L):
                u=L[a_]; v=L[m_]; pp=int((u&v).sum()); return pp,(pp*np_/(int(u.sum())*int(v.sum())) if pp else 0)
            pp,Rp=R(P); rd=[R(d) for d in D]; mean=np.mean([t[0] for t in rd])
            if pp<20 or mean<20: continue
            Rd=np.mean([t[1] for t in rd]); z=(Rp-Rd)/(Rd*math.sqrt(1/pp+1/(8*mean))); zz.append(z); out.append((dx,dy,pp,round(Rp,4),round(Rd,4),round(z,2)))
    print("%s: offsets tested %d, rms z %.2f, max |z| %.2f"%(label,len(zz),math.sqrt(np.mean(np.square(zz))),max(abs(v) for v in zz)))
    return out
A=run(0,N//2,"n <= 1e7"); B=run(N//2,N,"1e7 < n <= 2e7")
print("offset (dx,dy): prime pairs, ratio primes, ratio decoys, z   [first half | second half]")
for a,b in zip(A,B): print("  (%+d,%d): %7d %.4f %.4f %+5.2f | %7d %.4f %.4f %+5.2f"%(a[0],a[1],a[2],a[3],a[4],a[5],b[2],b[3],b[4],b[5]))
