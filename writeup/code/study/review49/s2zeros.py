# own number variance from the stored zeros (deduplicated), own unfolding
import numpy as np, json, math
from scipy.special import loggamma
Z=np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_zeros.npz')
FS=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']
LS=np.arange(0.5,10.01,0.5)
def theta(N,t):
    t=np.asarray(t,float); return 0.5*t*math.log(N)-t*math.log(2*math.pi)+np.imag(loggamma(1+1j*t))
def s2(u, ls=LS, step=0.01, lo=None, hi=None):
    # windows [x, x+l) with x on a grid in [u0, uend-l]
    out=[]
    a=u[0] if lo is None else lo; b=u[-1] if hi is None else hi
    for l in ls:
        x=np.arange(a, b-l, step)
        c=np.searchsorted(u,x+l)-np.searchsorted(u,x)
        out.append(c.var())
    return np.array(out)
res={}
for c in FS:
    g=np.unique(np.sort(Z[c['id']])); ndup=len(Z[c['id']])-len(g)
    for band,tag in (((20,400),'full'),((20,200),'lo'),((200,400),'hi'),((60,400),'b60')):
        b=g[(g>=band[0])&(g<=band[1])]; u=theta(c['N'],b)/math.pi
        res.setdefault(c['id'],{})[tag]=s2(u).tolist()
        res[c['id']]['n_'+tag]=len(b)
    res[c['id']]['ndup']=int(ndup)
json.dump(res,open('s2_own.json','w'))
R=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_readings.json'))
keep=[x['id'] for x in R['completeness'] if not x['excluded']]
A=np.array(R['sigma2']); B=np.array([res[k]['full'] for k in keep])
print('max |own - repo| over 136 curves x 20 l:', np.abs(A-B).max())
