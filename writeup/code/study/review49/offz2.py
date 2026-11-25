import json, glob, numpy as np, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
E={}
for f in glob.glob(A.D+'e2all/*.json'):
    r=json.load(open(f)); E[r['id']]=r
ids=sorted(i for i in E if i not in A.EXCL)
S=A.meas(ids); T=A.arr(ids,'full','full')
Tf=np.array([A.P[i]['full'] for i in ids]); Of=np.array([A.P[i]['offzeta'] for i in ids]); tm=Tf-Of
w=1/(S-Tf).var(0,ddof=1); rng=np.random.default_rng(3); n=len(ids)
def sl(y,x,ix=None):
    if ix is None: ix=np.arange(n)
    return A.slope(y[ix],x[ix],w)
prev=np.zeros_like(S)
for X in ('1000','10000','100000'):
    sX=np.array([E[i]['X'][X]['s2X'] for i in ids]); DX=A.arr(ids,'full','D'+X); c=sX-DX
    b=[sl(c,tm,rng.integers(0,n,n)) for _ in range(1000)]
    print('X=%-6s slope of [var(dS_X)-D_X] on offzeta template: %.3f +- %.3f'%(X,sl(c,tm),np.std(b)))
    # band between previous X and this X
    sR=np.array([E[i]['X'][X]['s2R'] for i in ids]); q=sR-(T-DX)
    b=[sl(q,tm,rng.integers(0,n,n)) for _ in range(1000)]
    print('         slope of remainder residual [var(dR_X)-(law-D_X)] on template: %.3f +- %.3f'%(sl(q,tm),np.std(b)))
b=[sl(S-Tf,tm,rng.integers(0,n,n)) for _ in range(1000)]
print('data - law on template (= gamma_offzeta - 1): %.3f +- %.3f'%(sl(S-Tf,tm),np.std(b)))
