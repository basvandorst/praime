import json, numpy as np, math, cypari2
from scipy.stats import spearmanr, pearsonr
pari=cypari2.Pari()
R=json.load(open('f5mine.json'))
mine={r['f']:r for r in R}
T={r['f']:r for r in json.load(open('/home/bas/praime/writeup/code/study/lowzero/r130.json'))['rows']}
print('max |g1 diff|',max(abs(mine[f]['g1']-T[f]['gamma1']) for f in T),'max |S diff|',max(abs(mine[f]['S']-T[f]['S']) for f in T))
# duplicates
from collections import defaultdict
byN=defaultdict(list)
for r in R: byN[r['N']].append(r['f'])
for N,fs in byN.items():
    for i in range(len(fs)):
        for j in range(i+1,len(fs)):
            if pari.nfisisom(pari(fs[i]),pari(fs[j])): print('ISOMORPHIC',N,fs[i],fs[j])
lnN=np.log([r['N'] for r in R]); S=np.array([r['S'] for r in R]); g1=np.array([r['g1'] for r in R])
g=g1*lnN/(2*np.pi)
def pred(w,loo=False,useS=True):
    out=[]
    for i in range(len(R)):
        m=np.abs(lnN-lnN[i])<w
        if loo: m[i]=False
        Sm=np.median(S[m]); gm=np.median(g1[m])
        rhs=(S[i]-Sm if useS else 0)+2/(0.25+gm**2)
        out.append(math.sqrt(max(2/rhs-0.25,0)) if rhs>0 else 10.0)
    return np.array(out)
for w in [0.3,0.5,1.0,2.0]:
    p=pred(w); pl=pred(w,True); p0=pred(w,useS=False)
    print('w',w,'spearman pred %.3f'%spearmanr(p,g1)[0],'LOO %.3f'%spearmanr(pl,g1)[0],'window-median-only %.3f'%spearmanr(p0,g1)[0],
          'pred-scaled vs g %.3f'%spearmanr(p*lnN,g1*lnN)[0], 'n rhs<=0',int(np.sum(p==10.0)))
print('spearman S vs g1',spearmanr(S,g1)[0],' lnN vs g1',spearmanr(lnN,g1)[0])
L1=np.array([r['L1'] for r in R]); print('spearman lnL1 vs g',spearmanr(np.log(L1),g)[0])
# residual: within-window S deviation vs gamma deviation
from numpy.polynomial import polynomial as Pn
S_res=S-np.polyval(np.polyfit(lnN,S,1),lnN); g_res=g1-np.polyval(np.polyfit(lnN,g1,1),lnN)
print('detrended spearman S vs g1',spearmanr(S_res,g_res)[0])
# fraction of S from first zero pair
c1=2/(0.25+g1**2); print('first pair share of S: median %.3f min %.3f'%(np.median(c1/S),np.min(c1/S)))
print('signatures',{str(k):sum(1 for r in R if r['mu']==k) for k in [[0,0,0,0],[0,0,1,1]]})
