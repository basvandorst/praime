import numpy as np, math, sys
from lib import sieve
N=20_000_000; s=sieve(N); y0=871; sv=np.ones(N+1,dtype=bool); sv[:2]=False
pr=np.nonzero(s[:int(N**0.5)+2])[0]
for q in pr[pr<=y0]: sv[2*q::q]=False
big=np.nonzero(s[:N//y0+2])[0]; big=big[big>y0]
RATIOS=((1,1),(2,1),(3,1),(3,2),(4,1),(5,1),(5,2),(5,3),(4,3),(5,4),(6,1),(6,5))
v2=sv.copy(); v3=sv.copy(); n2=n3=0; tot=0
for q in pr[pr>y0]:
    c=big[(big>=q)&(big<=N//q)].astype(np.int64); m=q*c; tot+=len(m)
    for (fa,fb) in RATIOS:
        U=fa*int(q)+fb*c; t=4*fa*fb*m; r=np.floor(np.sqrt(t.astype(np.float64))).astype(np.int64); r-=(r*r>t); r+=((r+1)*(r+1)<=t); hit=(U-r)<=2
        if (fa,fb)==(1,1): v2[m[hit]]=False
        v3[m[hit]]=False
print("composite survivors %d; struck with ratio 1: %d; with all ratios a/b, a <= 6: %d"%(tot,int(sv.sum()-v2.sum()),int(sv.sum()-v3.sum())))
n=np.arange(10**6,N-80000,dtype=np.int64)
def fl(x): return np.floor(x).astype(np.int64)
r=fl(np.sqrt(n)); r-=(r*r>n); r+=((r+1)**2<=n)
rules={"2 sqrt n":2*r,"4 sqrt n":4*r,"2 sqrt(2n)":2*fl(np.sqrt(2*n)),"sqrt(2n)":2*fl(np.sqrt(n/2)),"2 sqrt(3n)":2*fl(np.sqrt(3*n)),"2 sqrt(5n)":2*fl(np.sqrt(5*n)),"2 sqrt(6n)":2*fl(np.sqrt(6*n)),"3 sqrt n":2*fl(1.5*np.sqrt(n))}
print("root mean square deviation (%) of the pair ratio over k = -3..3:   sieve only | ratio 1 struck | ratios up to 6 struck | primes")
for name,h in rules.items():
    out=[]
    for L in (sv,v2,v3,s):
        a=int(L[n].sum()); v=[]
        for k in range(-3,4):
            m=n+h+2*k; v.append(int((L[n]&L[m]).sum())*len(n)/(a*int(L[m].sum())))
        mu=sum(v)/7; out.append(math.sqrt(sum((100*(x/mu-1))**2 for x in v)/7))
    print("  distance %-12s %.2f | %.2f | %.2f | %.2f"%(name,*out),flush=True)
