import json, glob, numpy as np, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
E={}
for f in glob.glob(A.D+'e2all/*.json'):
    r=json.load(open(f)); E[r['id']]=r
for lab,ids in (('136',sorted(i for i in E if i not in A.EXCL)),('138',sorted(E))):
    S=A.meas(ids); T=A.arr(ids,'full','full')
    Tf=np.array([A.P[i]['full'] for i in ids]); Of=np.array([A.P[i]['offzeta'] for i in ids])
    print(lab,'registered-style gamma_offzeta (frozen):  %.4f +- %.4f'%A.gam(S,Tf,Of,1000))
    print(lab,'own law, frozen template:                  %.4f +- %.4f'%A.gam(S,T,T-(Tf-Of),1000))
    for X in ('1000','10000','100000'):
        sX=np.array([E[i]['X'][X]['s2X'] for i in ids]); DX=A.arr(ids,'full','D'+X)
        Tc=T-DX+sX
        print(lab,'corrected law X=%-6s gamma_offzeta %.4f +- %.4f   gamma_k1p3 %.4f +- %.4f  gamma_w1p3 %.4f +- %.4f'%(X,*A.gam(S,Tc,Tc-(Tf-Of),1000),*A.gam(S,Tc,Tc+A.arr(ids,'full','dk1p3'),1000),*A.gam(S,Tc,Tc+A.arr(ids,'full','dw1p3'),1000)))
    for X in ('1000','H'):
        EF=A.arr(ids,'full','D'+X)+A.arr(ids,'full','G'+X)
        w=1/(S-Tf).var(0,ddof=1)
        print(lab,'pseudo-data EF_%s gamma_offzeta %.3f'%(X,A.slope(EF-Of,Tf-Of,w)))
