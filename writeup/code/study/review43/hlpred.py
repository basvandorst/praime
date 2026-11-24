# Exact HL-model coefficient for E[chi2]/dof = 1 + a(q,H)/lam  (finite-H sums), and asymptotic law a_asym.
import numpy as np, math, json
from sympy import factorint, totient
HMAX=10**7
C2=0.66016181584686957392781211
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
P=primes_upto(HMAX)
G=np.ones(HMAX)
for p in P[1:]:
    G[p::p]*=(p-1)/(p-2)
S2=np.zeros(HMAX); S2[0::2]=2.0
gamma=0.57721566490153286; B=1-gamma-math.log(2*math.pi)
e=np.arange(HMAX)
def Sprime(q):
    f=factorint(q); Sp=C2*G.copy()
    if 2 not in f: Sp*=S2
    for p in f:
        if p==2: continue
        div=(e%p==0); Sp[div]/=p/(p-1); Sp[~div]/=p*(p-2)/(p-1)**2
    return Sp
out={}
for q in [3,4,5,7,8,9,10,12,30,210]:
    Sp=Sprime(q); phi=int(totient(q)); f=factorint(q)
    Sq=sum(math.log(p)/(p-1) for p in f)
    red=[a for a in range(q) if math.gcd(a,q)==1]
    for H in [10**4,10**5,10**6,10**7]:
        T=np.zeros(q)
        for c in range(q):
            idx=np.arange(c if c>0 else q, H, q)
            T[c]=np.sum((H-idx)*(Sp[idx]-1.0))
        A=np.zeros((phi,phi))
        for i,a in enumerate(red):
            for j,b in enumerate(red):
                c=(b-a)%q
                A[i,j]=(q/(phi*H))*(T[c]+T[(-c)%q]) - (q/phi if i==j else 0)
        a_exact=(np.trace(A)-A.sum()/phi)/(phi-1)
        a_asym=-(math.log(H/q)+Sq-B) + (math.log(q)-Sq)/(phi-1)
        vq_coef=A[0,0]  # Var/mean = 1 + A00/lam
        out['%d_%d'%(q,H)]=dict(a_exact=a_exact,a_asym=a_asym,vq_exact=A[0,0],vq_asym=-(math.log(H/q)+Sq-B),A=A.tolist())
        print(q,H,'a_exact %.4f a_asym %.4f  vq %.4f vs %.4f'%(a_exact,a_asym,A[0,0],-(math.log(H/q)+Sq-B)))
    # q=1 total
for H in [10**4,10**5,10**6,10**7]:
    h=np.arange(1,H); Sfull=S2[1:H]*C2*G[1:H]
    T=np.sum((H-h)*(Sfull-1)); out['1_%d'%H]=dict(v1=2*T/H-1, v1_asym=-(math.log(H)-B))
json.dump(out,open('hlpred.json','w'))
