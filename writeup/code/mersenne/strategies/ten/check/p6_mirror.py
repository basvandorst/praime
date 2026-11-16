import numpy as np
from lib import sieve
N=10**7; s=sieve(N); 
def rev_arr(n,B,d):
    r=np.zeros_like(n); x=n.copy()
    for _ in range(d): r=r*B+x%B; x//=B
    return r
def study(B,d,qs):
    lo,hi=B**(d-1),min(B**d,N+1); P=np.nonzero(s[lo:hi])[0]+lo; r=rev_arr(P,B,d)
    ok=np.gcd(r,B)==1; P=P[ok]; r=r[ok]
    K2=B/len([a for a in range(B) if np.gcd(a,B)==1])
    for q in range(2,B+2):
        if all(q%t for t in range(2,q)) and (B*B-1)%q==0 and B%q: K2*=q/(q-1)
    base=float(np.sum(K2/np.log(r[r>1]))); corr=1.0; out=[]
    for q in qs:
        m=(P%q!=0)&(r!=q); frac=np.mean(r[m]%q==0); c=(1-frac)/(1-1/q); corr*=c; out.append("q=%d: %.4f vs %.4f"%(q,frac,1/q))
    if B**d<=10**7*B and r.max()<=N*0+10**9:
        from sympy import isprime
        real=sum(1 for x in r.tolist() if (s[x] if x<=N else isprime(x)))
    print("base %d d=%d primes(rev coprime) %d  real %d  plain %.0f  with q-corrections %.0f  (factor %.4f)"%(B,d,len(P),real,base,base*corr,corr)); print("   P(q | rev) measured vs 1/q:", "; ".join(out))
study(10,6,[7,13,37,41,101,271,17,19]); study(10,7,[7,13,37,41,101,271,17,19]); study(10,5,[7,13,37,41,101,271,17,19])
study(2,15,[5,7,17,31,73,127,11,13]); study(2,18,[5,7,17,31,73,127,11,13]); study(2,20,[5,7,17,31,73,127,11,13])
