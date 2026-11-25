import json, glob, numpy as np, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
E={}
for f in glob.glob(A.D+'e2all/*.json'):
    r=json.load(open(f)); E[r['id']]=r
ids=sorted(i for i in E if i not in A.EXCL)
S=A.meas(ids); T=A.arr(ids,'full','full'); r=S-T; pick=[0,1,3,7,11,15,19]
for X in ('1000','10000','100000'):
    sX=np.array([E[i]['X'][X]['s2X'] for i in ids]); DX=A.arr(ids,'full','D'+X); c=sX-DX
    cc=np.array([np.corrcoef(r[:,j],c[:,j])[0,1] for j in range(20)])
    fr=1-((r-c).var(0)/r.var(0))
    print('X=%-6s corr(data-law, var(dS_X)-D_X) %s ; share of residual variance removed %s'%(X,A.fmt(cc[pick],'%+.2f'),A.fmt(fr[pick],'%+.2f')))
