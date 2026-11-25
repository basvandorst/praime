import os, numpy as np, json, math
R={}; R.update(json.load(open('remap_B.json'))); R.update(json.load(open('remap_Dn_Fn_C_Bp_R120.json')))
P=json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/spacing/predictions_spacing.json'))
PR=(2,3,5,7,11,13); EXCL={'Fn':[1153],'C':[197,-395]}
def design(D): return np.array([[1.0]+[float(abs(d)%p==0) for p in PR] for d in D])
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
summary={}
for name in ['B','Dn','Fn','C','Bp','R120']:
    S=P[name]; Ds=[d for d in S['D'] if d not in EXCL.get(name,[])]; idx={d:i for i,d in enumerate(S['D'])}; ii=[idx[d] for d in Ds]
    X=design(Ds); n=len(Ds); rec=[R[name][str(d)] for d in Ds]
    rng=np.random.default_rng(3); BI=[rng.integers(0,n,n) for _ in range(1000)]
    summary[name]={}
    for k,stn in ((0,'var'),(1,'r')):
        y=np.array([r['v0' if k==0 else 'r0'] for r in rec])
        ign=np.array(S['pred']['ignored'][stn])[ii]; law=np.array(S['pred']['law'][stn])[ii]
        c=ols(X,y-ign); cl=ols(X,law-ign)
        line=[]
        for j,p in enumerate(PR):
            A=np.array([abs(d)%p==0 for d in Ds])
            add=np.array([rec[i]['add%d'%p][k] if A[i] else 0 for i in range(n)])
            adds=np.array([rec[i]['adds%d'%p][k] if A[i] else 0 for i in range(n)])
            remt=np.array([rec[i]['remt%d'%p][k] if not A[i] else 0 for i in range(n)])
            remr=np.array([rec[i]['remr%d'%p][k] if not A[i] else 0 for i in range(n)])
            # regression-matched predicted coefficient from the real base: c - c_after
            pa=c[j+1]-ols(X,y-ign+add)[j+1]      # ADD: surrogate replaces p|D
            pt=ols(X,y-ign+remt)[j+1]-c[j+1]     # REMT: true-wave removal on p not| D -> coefficient change = -(predicted c)
            pt=-pt
            pas=c[j+1]-ols(X,y-ign+adds)[j+1]
            pr=-(ols(X,y-ign+remr)[j+1]-c[j+1])
            # bootstrap se of c - predicted (ADD and REMT)
            dA=np.array([ (ols(X[b],(y-ign)[b])[j+1]) - ( ols(X[b],(y-ign)[b])[j+1]-ols(X[b],(y-ign+add)[b])[j+1]) for b in BI])
            dT=np.array([ ols(X[b],(y-ign+remt)[b])[j+1] for b in BI])
            line.append(dict(p=p,c=c[j+1],law=cl[j+1],add=pa,adds=pas,remt=pt,remr=pr,se_after_add=dA.std(),se_after_remt=dT.std(),
                             add_after=c[j+1]-pa, remt_after=c[j+1]-pt,
                             add_raw=-add[A].mean(),remt_raw=remt[~A].mean(),remr_raw=remr[~A].mean(),adds_raw=-adds[A].mean()))
        summary[name][stn]=line
        print('%-4s %-3s'%(name,stn))
        print('    p    data c     law   realADD  realREMT  randREM  ADDflip | ADD/law REMT/law | data-ADD (se)  data-REMT (se)')
        for L in line:
            print('   %2d %+.5f %+.5f %+.5f %+.5f %+.5f %+.5f | %6.3f %6.3f | %+.5f (%.5f) %+.5f (%.5f)'%(L['p'],L['c'],L['law'],L['add'],L['remt'],L['remr'],L['adds'],L['add']/L['law'],L['remt']/L['law'],L['add_after'],L['se_after_add'],L['remt_after'],L['se_after_remt']))
json.dump(summary,open('remap_summary.json','w'))
