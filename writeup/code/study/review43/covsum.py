# Evaluate HL restricted sums directly. Own code.
import numpy as np, math
from sympy import factorint, totient
HMAX=10**7
S=np.load('S_%d.npy'%HMAX)
C2=0.66016181584686957392781211
G=np.ones(HMAX); 
# G(e)=prod_{p odd | e}(p-1)/(p-2)
S2=np.zeros(HMAX); S2[0::2]=2.0
Gs = S.copy(); Gs[0::2] = S[0::2]/(2*C2)
# odd e: need G separately; recompute G for all
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
P=primes_upto(HMAX)
for p in P[1:]:
    G[p::p]*=(p-1)/(p-2)
gamma=0.57721566490153286; B=1-gamma-math.log(2*math.pi)
def Sprime(q):
    f=factorint(q)
    Sp = C2*G.copy()
    if 2 not in f: Sp*=S2
    e=np.arange(HMAX)
    for p in f:
        if p==2: continue
        div = (e%p==0)
        Sp[div] /= p/(p-1)
        Sp[~div] /= p*(p-2)/(p-1)**2
    return Sp
def Lam(m):
    f=factorint(m)
    if len(f)==1: return math.log(list(f)[0])
    return 0.0
print('B =',B)
# q=1 check
for H in [10**4,10**5,10**6,10**7]:
    h=np.arange(1,H)
    T=np.sum((H-h)*(S[1:H]-1.0))
    print('q=1 H=%d  2T/H + lnH = %.5f'%(H, 2*T/H+math.log(H)))
for q in [3,4,5,7,8,9,10,12,30,210]:
    Sp=Sprime(q); phi=int(totient(q)); f=factorint(q)
    Sq=sum(math.log(p)/(p-1) for p in f)
    for H in [10**4,10**5,10**6,10**7]:
        e=np.arange(H); w=(H-e).astype(float)
        T={}
        for c in range(q):
            idx=np.arange(c if c>0 else q, H, q)
            T[c]=np.sum((H-idx)*(Sp[idx]-1.0))
        D=2*q/(phi*H)*T[0]
        # D should equal q/phi - ln(H/q) - Sq + B if law includes diagonal p^2 term
        line='q=%d H=%.0e  D+ln(H/q)+Sq = %.4f  (q/phi+B=%.4f, B=%.4f) |'%(q,H,D+math.log(H/q)+Sq, q/phi+B, B)
        offs=[]
        for c in range(1,q):
            # only c that connect reduced classes: exists a reduced with a+c reduced
            ok=any(math.gcd(a,q)==1 and math.gcd(a+c,q)==1 for a in range(q))
            if not ok: continue
            A=q/(phi*H)*(T[c]+T[(-c)%q])
            m=q//math.gcd(c,q)
            offs.append((c, round(A,4), round(-Lam(m)/int(totient(m)),4)))
        if H in (10**5,10**7) or q in (4,9,210):
            print(line, offs if q<=12 else [(c,a,b) for c,a,b in offs if abs(a-b)>0.003][:8], 'maxdev', max(abs(a-b) for c,a,b in offs))
        else:
            print(line, 'maxdev', max(abs(a-b) for c,a,b in offs))
