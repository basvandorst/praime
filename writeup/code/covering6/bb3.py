# Branch-and-bound for the level-n tail relaxation: hierarchical branching + symmetry; writes an exact certificate tree.
# Node: for each base modulus m a restriction  r_m = rho_m (mod g_m),  g_m | m.   m is "fixed" when g_m = m.
# U = residues mod L outside the fixed classes.
# Leaf {"w":{b:W_b}}: integers W_b >= 0 with
#   N := W(U) - sum_{m not fixed} max_{r = rho_m mod g_m} W(U & class r mod m)  >=  D := sum_{d|45} max_r W(U & class r mod 2^n d) > 0.
# Internal node {"m":m,"q":q,"ch":[...]} (q prime, g_m*q | m): child j has g_m -> g_m*q, rho_m -> rho_m + g_m*j.
# A child may be {"sym":[j2,perm]}: perm is a permutation of Z/q^e (the q-component of Z/L under CRT) that preserves
# congruence modulo every power of q, maps every restriction set of the node onto itself and child j onto child j2.
import sys, time, json, itertools
import numpy as np
from fractions import Fraction
from math import lcm
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
DS=(1,3,5,9,15,45)
n=int(sys.argv[1]); TL=float(sys.argv[2]); out=sys.argv[3]; ORDER=sys.argv[4] if len(sys.argv)>4 else 'size'
L=2**n*45
M=[m for m in range(6,L+1) if L%m==0]
TOP=[2**n*d for d in DS]
PP={2:2**n,3:9,5:5}
def ppart(m,p):
    k=1
    while m%(k*p)==0 and PP[p]%(k*p)==0: k*=p
    return k
def gen_group(p):
    N=PP[p]; divs=[]; k=p
    while k<N: divs.append(k); k*=p
    G=[]
    for pi in itertools.permutations(range(N)):
        ok=True
        for dv in divs:
            img={}
            for x in range(N):
                if img.setdefault(x%dv,pi[x]%dv)!=pi[x]%dv: ok=False;break
            if not ok: break
        if ok: G.append(pi)
    return G
GRP={p:gen_group(p) for p in (2,3,5)}
print('group sizes',{p:len(GRP[p]) for p in GRP},flush=True)
def compset(p,g,rho):
    gp=ppart(g,p)
    return frozenset(x for x in range(PP[p]) if x%gp==rho%gp)
if ORDER=='size': RANK={m:i for i,m in enumerate(M)}
else:
    RANK={m:i for i,m in enumerate(sorted(M,key=lambda m:((m&-m),m)))}
stats={'lp':0,'leaf':0,'int':0,'sym':0,'maxdepth':0,'fail_round':0,'frac':0.0}
t0=time.time()
def uncovered(st):
    cov=np.zeros(L,bool)
    for m,(g,rho) in st.items():
        if g==m: cov[rho::m]=True
    return np.where(~cov)[0]
def exact_check(st,W):
    U=[int(b) for b in uncovered(st)]
    N=sum(W.get(b,0) for b in U)
    for m,(g,rho) in st.items():
        if g==m: continue
        s=[0]*m
        for b in U: s[b%m]+=W.get(b,0)
        N-=max(s[r] for r in range(rho,m,g))
    Dn=0
    for m in TOP:
        s=[0]*m
        for b in U: s[b%m]+=W.get(b,0)
        Dn+=max(s)
    return N,Dn
def lp(st):
    U=uncovered(st); nu=len(U)
    free=[m for m in M if st[m][0]!=m]
    nv=nu+len(free)+6
    rows=[];cols=[];vals=[]; nr=0
    for j,m in enumerate(free):
        g,rho=st[m]; res=U%m
        for r in range(rho,m,g):
            ii=np.where(res==r)[0]
            if len(ii)==0: continue
            rows+= [nr]*(len(ii)+1); cols+=list(ii)+[nu+j]; vals+=[1.0]*len(ii)+[-1.0]; nr+=1
    for j,m in enumerate(TOP):
        res=U%m
        for r in np.unique(res):
            ii=np.where(res==r)[0]
            rows+= [nr]*(len(ii)+1); cols+=list(ii)+[nu+len(free)+j]; vals+=[1.0]*len(ii)+[-1.0]; nr+=1
    A=csr_matrix((vals,(rows,cols)),shape=(nr,nv))
    c=np.zeros(nv); c[:nu]=-1; c[nu:nu+len(free)]=1
    Aeq=np.zeros((1,nv)); Aeq[0,nu+len(free):]=1
    r=linprog(c,A_ub=A,b_ub=np.zeros(nr),A_eq=Aeq,b_eq=[1],bounds=(0,None),method='highs-ds')
    stats['lp']+=1
    return -r.fun, U, r.x[:nu], free
def rationalise(st,U,w):
    if len(w)==0 or max(w)<=0: return None
    for den in (1,2,3,4,5,6,8,9,10,12,15,18,20,24,30,36,40,45,60,72,90,120,180,360,720,2520,5040):
        fr=[Fraction(float(v)).limit_denominator(den*50) for v in w]
        l=1
        for f in fr: l=lcm(l,f.denominator)
        if l>10**9: continue
        W={int(b):int(f*l) for b,f in zip(U,fr) if f>0}
        if not W: continue
        N,Dn=exact_check(st,W)
        if Dn>0 and N>=Dn: return W
    return None
def rec(st,stab,depth,vol=1.0):
    stats['maxdepth']=max(stats['maxdepth'],depth)
    if time.time()-t0>TL: raise TimeoutError
    v,U,w,free=lp(st)
    if stats['lp']%500==0: print('  ..',stats,'time %.0f'%(time.time()-t0),flush=True)
    if v>=1-1e-7:
        W=rationalise(st,U,w)
        if W is not None:
            stats['leaf']+=1; stats['frac']+=vol
            return {'w':W}
        stats['fail_round']+=1
    if not free:
        print('COUNTEREXAMPLE or rounding failure at full assignment',st,'LP value',v,flush=True)
        return {'bad':{m:st[m][1] for m in st},'v':v}
    m=min(free,key=lambda m:RANK[m]); g,rho=st[m]
    q=[p for p in (2,3,5) if (m//g)%p==0][0]
    stats['int']+=1
    node={'m':m,'q':q,'ch':[None]*q}
    csets=[compset(q,g*q,rho+g*j) for j in range(q)]
    idx={s:j for j,s in enumerate(csets)}
    rep=list(range(q)); via={}
    for pi in stab[q]:
        for j in range(q):
            j2=idx[frozenset(pi[x] for x in csets[j])]
            if j2<rep[j]: rep[j]=j2; via[j]=pi
    # make representatives consistent (minimum of orbit)
    for j in range(q):
        if rep[j]!=j and rep[rep[j]]!=rep[j]:
            raise RuntimeError('orbit rep inconsistency')
    nrep={j:sum(1 for i in range(q) if rep[i]==j) for j in range(q) if rep[j]==j}
    for j in range(q):
        if rep[j]==j:
            s2=dict(st); s2[m]=(g*q,rho+g*j)
            stab2=dict(stab); stab2[q]=[pi for pi in stab[q] if frozenset(pi[x] for x in csets[j])==csets[j]]
            node['ch'][j]=rec(s2,stab2,depth+1,vol*nrep[j]/q)
        else:
            node['ch'][j]={'sym':[rep[j],list(via[j])]}; stats['sym']+=1
    if depth<=4: print('depth',depth,'branch',m,q,'orbits',sorted(nrep.items()),'done; stats',stats,'time %.0f'%(time.time()-t0),flush=True)
    return node
try:
    tree=rec({m:(1,0) for m in M},dict(GRP),0)
    json.dump({'n':n,'tree':tree},open(out,'w'))
    print('FINISHED',stats,'time %.0f'%(time.time()-t0),flush=True)
except TimeoutError:
    print('TIMEOUT',stats,flush=True)
