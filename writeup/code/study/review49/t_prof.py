import sys, json, time, numpy as np, math
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mylaw as M
f=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves'][0]
t0=time.time(); E=M.pari.ellinit(f['ainv']); an=M.pari.ellan(E,M.PMAX); a=np.array([int(x) for x in an]); print('ellan', time.time()-t0)
t0=time.time(); Ls=M.pari.lfunsympow(E,2); Li=M.pari.lfuninit(Ls,[2,0,8.3],2); print('lfuninit', time.time()-t0)
t0=time.time()
for e in M.EPS[:200]:
    s=M.pari('2+%.10f*I'%e); M.pari.lfun(Li,s,0); M.pari.lfun(Li,s,1); M.pari.lfun(Li,s,2)
print('200 eps x3', time.time()-t0)
