# independent: (A) order-2 Bahadur per class mod 30030; (B) exact HL inclusion-exclusion with full singular series
import math, sys, itertools, json
import numpy as np
from sympy import primerange
lam = math.log(1.5e9)
SMALL=[2,3,5,7,11,13]; PP=30030
PR = np.array(list(primerange(17, 2_000_000)), dtype=np.float64)
QX = 2000  # exact residue counts up to QX, all points distinct beyond (spans < 1300)
EX = [int(q) for q in PR if q <= QX]; GEN = PR[PR > QX]
# tail beyond 2e6: sum_{q>Q} [ln(1-k/q) - k ln(1-1/q)] ~ -k(k-1)/2 * sum 1/q^2 ~ -k(k-1)/2/(Q ln Q)
Qe = 2e6
GT = {}
def gtail(k):
    if k not in GT:
        GT[k] = float(np.sum(np.log1p(-k/GEN) - k*np.log1p(-1/GEN))) - k*(k-1)/2/(Qe*math.log(Qe))
    return GT[k]
def lnS(pts, qs):
    """ln of prod over q in qs (exact) and q > QX (generic) of (1-nu/q)/(1-1/q)^k"""
    k = len(pts); s = 0.0
    for q in qs:
        nu = len({x % q for x in pts})
        if nu >= q: return -math.inf
        s += math.log(1 - nu/q) - k*math.log(1 - 1/q)
    return s + gtail(k)
def lnS13(pts): return lnS(pts, EX)
def lnSfull(pts): return lnS(pts, [2,3,5,7,11,13] + EX)

def order2(d, order=2):
    B = [0,2,d,d+2]; H = list(range(6, d, 6)); n=len(H)
    lB = lnS13(B); l02 = lnS13([0,2])
    ls = np.array([lnS13(B+[h,h+2]) - lB for h in H])
    R = np.zeros((n,n))
    if order == 2:
        for i in range(n):
            for j in range(i+1,n):
                R[i,j]=R[j,i]=math.expm1(lnS13(B+[H[i],H[i]+2,H[j],H[j]+2]) + lB - (ls[i]+lB) - (ls[j]+lB))
    rho0 = 1/(lam*math.prod(1-1/q for q in SMALL))
    p = rho0**2*np.exp(ls); w = p/(1-p)
    tot = 0.0; ncls = 0
    for a in range(PP):
        if math.gcd(a,PP)!=1 or math.gcd(a+2,PP)!=1: continue
        ncls += 1
        if math.gcd(a+d,PP)!=1 or math.gcd(a+d+2,PP)!=1: continue
        idx = [i for i,h in enumerate(H) if math.gcd(a+h,PP)==1 and math.gcd(a+h+2,PP)==1]
        L = float(np.sum(np.log1p(-p[idx])))
        if order==2 and len(idx)>1:
            ii=np.array(idx); L += 0.5*float(w[ii]@R[np.ix_(ii,ii)]@w[ii])
        tot += rho0**2*math.exp(lB-l02)*math.exp(L)
    return tot/ncls, ncls

def exact(d):
    B=[0,2,d,d+2]; H=list(range(6,d,6)); l02=lnSfull([0,2])
    tot=0.0
    for m in range(len(H)+1):
        for T in itertools.combinations(H,m):
            pts = B+[x for h in T for x in (h,h+2)]
            l = lnSfull(pts)
            if l == -math.inf: continue
            tot += (-1)**m * math.exp(l - l02 - (2+2*m)*math.log(lam))
    return tot

if __name__=='__main__':
    ref = json.load(open('/home/bas/praime/writeup/data/study/r59.json'))['law']
    for d in [int(x) for x in sys.argv[1].split(',')]:
        o2,nc = order2(d); o1,_ = order2(d,1)
        line = 'd=%d ncls=%d order1=%.8f order2=%.8f repo=%.8f rel(o2/repo-1)=%.2e' % (d,nc,o1,o2,ref[str(d)][1],o2/ref[str(d)][1]-1)
        if d <= int(sys.argv[2]):
            ex = exact(d); line += ' exact=%.8f o2/exact-1=%.2e o1/exact-1=%.2e' % (ex, o2/ex-1, o1/ex-1)
        print(line, flush=True)
