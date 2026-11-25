import os, numpy as np, json, math
from sympy import primefactors
st=json.load(open('base_stats.json'))
P=json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/spacing/predictions_spacing.json'))
PR=(2,3,5,7,11,13)
def design(D): return np.array([[1.0]+[float(abs(d)%p==0) for p in PR] for d in D])
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
for a,b in (('B','Bp'),('R120','C')):
    common=sorted(set(P[a]['D'])&set(P[b]['D'])); X=design(common)
    lq=np.log(np.abs(common)); lp=np.array([sum(1/p for p in primefactors(abs(d)) if p>13) for d in common])
    sg=np.array([1.0 if d>0 else 0.0 for d in common])
    for stn in ('var','r'):
        out=[]
        for nm in (a,b):
            S=P[nm]; idx={d:i for i,d in enumerate(S['D'])}; ii=[idx[d] for d in common]; rec={x['D']:x for x in st[nm]}
            y=np.array([rec[d][stn] for d in common]); law=np.array(S['pred']['law'][stn])[ii]; ign=np.array(S['pred']['ignored'][stn])[ii]
            raw=y-X@ols(X,y); res=(y-law)-X@ols(X,y-law); lawres=(law-ign)-X@ols(X,law-ign)
            out.append((raw,res,lawres))
        print(a,b,stn,'corr raw-resid %.2f | data-law resid %.2f | law-ign resid (law structure beyond indicators) %.2f'%(np.corrcoef(out[0][0],out[1][0])[0,1],np.corrcoef(out[0][1],out[1][1])[0,1],np.corrcoef(out[0][2],out[1][2])[0,1]))
        for k,nm in enumerate((a,b)):
            print('    %-4s corr(resid, ln|D|) raw %.2f d-l %.2f | corr(resid, large-prime load) raw %+.2f d-l %+.2f law %+.2f | corr(resid, D>0) raw %+.2f d-l %+.2f  | sd raw %.5f d-l %.5f lawres %.5f'%(
                nm,np.corrcoef(out[k][0],lq)[0,1],np.corrcoef(out[k][1],lq)[0,1],np.corrcoef(out[k][0],lp)[0,1],np.corrcoef(out[k][1],lp)[0,1],np.corrcoef(out[k][2],lp)[0,1],
                np.corrcoef(out[k][0],sg)[0,1],np.corrcoef(out[k][1],sg)[0,1],out[k][0].std(),out[k][1].std(),out[k][2].std()))
