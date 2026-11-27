# Tail relaxation at level n.
#   x_{m,r} in {0,1}: one class for each modulus m | 2^n*45, m >= 6   (the real congruences up to level n)
#   y_{d,r} >= 0 real: fractional classes mod 2^n*d, d | 45           (all levels > n, merged)
#   every b mod 2^n*45:  sum_m x_{m,b mod m} + sum_d y_{d,b mod 2^n d} >= 1
#   every d: sum_r y_{d,r} <= t.      Minimise t.
# A covering with LCM 2^a*45 (a > n) gives a feasible point with t = 1 - 2^(n-a) < 1.
# So min t >= 1 at one level n proves that no covering exists for any a.
import sys, time, json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix
DS = (1,3,5,9,15,45)
def base_moduli(n, mmin=6):
    return sorted(2**k*d for k in range(n+1) for d in DS if 2**k*d >= mmin)
def build(n, sym=True, mmin=6):
    L=2**n*45; M=base_moduli(n,mmin); idx={}
    for m in M:
        for r in range(m): idx[('x',m,r)]=len(idx)
    nx=len(idx)
    for d in DS:
        for r in range(2**n*d): idx[('y',d,r)]=len(idx)
    nv=len(idx)+1
    rows=[];cols=[]
    for b in range(L):
        for m in M: rows.append(b); cols.append(idx[('x',m,b%m)])
        for d in DS: rows.append(b); cols.append(idx[('y',d,b%(2**n*d))])
    A1=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(L,nv))
    rows=[];cols=[]
    for i,m in enumerate(M):
        for r in range(m): rows.append(i); cols.append(idx[('x',m,r)])
    A2=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(M),nv))
    rows=[];cols=[];vals=[]
    for i,d in enumerate(DS):
        for r in range(2**n*d): rows.append(i); cols.append(idx[('y',d,r)]); vals.append(1.0)
        rows.append(i); cols.append(nv-1); vals.append(-1.0)
    A3=csr_matrix((vals,(rows,cols)),shape=(6,nv))
    lb=np.zeros(nv); ub=np.ones(nv); ub[-1]=10
    if sym:
        lb[idx[('x',9,0)]]=1
        for r in range(15):
            if r%5: ub[idx[('x',15,r)]]=0
        for r in range(L):
            if r%2**n: ub[idx[('x',L,r)]]=0
    integ=np.zeros(nv); integ[:nx]=1
    c=np.zeros(nv); c[-1]=1
    return L,M,idx,nx,nv,c,[LinearConstraint(A1,1,np.inf),LinearConstraint(A2,-np.inf,1),LinearConstraint(A3,-np.inf,0)],integ,Bounds(lb,ub)
if __name__=='__main__':
    n=int(sys.argv[1]); tl=float(sys.argv[2]); t0=time.time()
    L,M,idx,nx,nv,c,cons,integ,bnds=build(n)
    res=milp(c,constraints=cons,integrality=integ,bounds=bnds,options={'time_limit':tl,'mip_rel_gap':0.0,'disp':len(sys.argv)>3})
    print('n',n,'L',L,'status',res.status,'min t',res.fun,'dual bound',getattr(res,'mip_dual_bound',None),'time %.1f'%(time.time()-t0))
    if res.x is not None:
        X={m:[r for r in range(m) if res.x[idx[('x',m,r)]]>0.5] for m in M}
        Y={d:{r:round(float(res.x[idx[('y',d,r)]]),6) for r in range(2**n*d) if res.x[idx[('y',d,r)]]>1e-9} for d in DS}
        print('X',json.dumps(X)); print('Y',json.dumps(Y))
