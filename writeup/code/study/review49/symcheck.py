import json, numpy as np, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
ids=sorted(i for i in A.LAW if i in A.S2 and i not in A.EXCL)
S=A.meas(ids); T=A.arr(ids,'full','full')
Tf=np.array([A.P[i]['full'] for i in ids]); Of=np.array([A.P[i]['offzeta'] for i in ids]); K1f=np.array([A.P[i]['k1'] for i in ids])
nos=A.arr(ids,'full','nosym')
w=1/(S-Tf).var(0,ddof=1)
def sl(y,x): return A.slope(y,x,w)
print('registered gamma_offzeta on data (frozen):', round(sl(S-Of,Tf-Of),4))
print('pseudo-data = own law      : gamma_offzeta', round(sl(T-Of,Tf-Of),4), ' gamma_nosym', round(sl(T-nos,T-nos),4))
for X in ('300','1000','10000','H'):
    EF=A.arr(ids,'full','D'+X)+A.arr(ids,'full','G'+X)
    print('pseudo-data = EF_%-5s     : gamma_offzeta %.3f   gamma_nosym %.3f   gamma_A %.3f'%(X, sl(EF-Of,Tf-Of), sl(EF-nos,T-nos), sl(EF,Tf)))
print('template sd (law - offzeta, frozen):', A.fmt((Tf-Of).std(0)[[1,3,7,11,15,19]],'%.4f'))
print('template sd (law - nosym, own)    :', A.fmt((T-nos).std(0)[[1,3,7,11,15,19]],'%.4f'))
print('corr per l between (law-EF_1000) and (law-offzeta):', A.fmt(np.array([np.corrcoef((T-A.arr(ids,'full','D1000')-A.arr(ids,'full','G1000'))[:,j],(Tf-Of)[:,j])[0,1] for j in range(20)])[[1,3,7,11,15,19]],'%+.2f'))
