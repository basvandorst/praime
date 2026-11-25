import json, numpy as np, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
E=json.load(open('e2_out.json'))
pick=[0,1,3,7,11,15,19]
ids=[e['id'] for e in E]
S=np.array([A.S2[i]['full'] for i in ids]); T=A.arr(ids,'full','full')
meas=np.array([e['meas'] for e in E]); print('meas(true S) - registered-estimator S2 max', np.abs(meas-S).max())
for X in ['10','100','1000','10000','100000']:
    sX=np.array([e['X'][X]['s2X'] for e in E]); sR=np.array([e['X'][X]['s2R'] for e in E])
    DX=A.arr(ids,'full','D'+X)
    cross=S-sX-sR
    print('X=%-6s var(dS_X) - D_X: mean %s  sd %s'%(X,A.fmt((sX-DX).mean(0)[pick]),A.fmt((sX-DX).std(0)[pick],'%.4f')))
    print('         data - var(dS_X): mean %s ; var(dR_X) %s ; 2cov %s'%(A.fmt((S-sX).mean(0)[pick],'%+.3f'),A.fmt(sR.mean(0)[pick],'%+.3f'),A.fmt(cross.mean(0)[pick],'%+.4f')))
    # does var(dS_X)-D_X track the per-curve residual data - law ?
    r=S-T; c=[np.corrcoef(r[:,j],(sX-DX)[:,j])[0,1] for j in range(20)]
    print('         corr over curves of (data - law) with (var(dS_X) - D_X):', A.fmt(np.array(c)[pick],'%+.2f'))
print()
for X in ['100','1000','10000','100000']:
    sX=np.array([e['X'][X]['s2X'] for e in E]); sR=np.array([e['X'][X]['s2R'] for e in E]); DX=A.arr(ids,'full','D'+X)
    q=sR-(T-DX); n=len(ids)
    print('X=%-6s remainder test: mean[var(dR_X) - (law - D_X)] %s'%(X,A.fmt(q.mean(0)[pick])))
    print('          se %s   (data - law on the same curves: se %s, mean %s)'%(A.fmt(q.std(0,ddof=1)[pick]/np.sqrt(n),'%.4f'),A.fmt((S-T).std(0,ddof=1)[pick]/np.sqrt(n),'%.4f'),A.fmt((S-T).mean(0)[pick])))
