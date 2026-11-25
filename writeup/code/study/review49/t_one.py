import sys, json, time, numpy as np, math
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mylaw as M
P=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_predictions.json'))
print('GUE own - frozen', np.abs(M.GUE-np.array(P['gue'])).max())
c=P['curves'][0]
t0=time.time(); cv=M.Curve(c['ainv'] if len(c['ainv'])==5 else None, c['N']) if False else None
FS={x['id']:x for x in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
f=FS[c['id']]
t0=time.time(); cv=M.Curve(f['ainv'], f['N']); print('init', time.time()-t0, 'LS1', cv.LS1, 'frozen', c['LS1'])
tn,wn=M.nodes(n=12); Ln=M.Lt(f['N'],tn); wt=wn*Ln; wt/=wt.sum()
FJ=cv.FJ(); FK=cv.FK()
law=np.zeros(20)
for L,w in zip(Ln,wt):
    law+=w*(M.GUE+M.s2_from_F(FJ,FK,L)+cv.diag_corr(L))
print('time', time.time()-t0)
print('own law', np.round(law,5)[[0,1,3,7,11,15,19]])
print('frozen ', np.round(np.array(c['full']),5)[[0,1,3,7,11,15,19]])
print('max diff', np.abs(law-np.array(c['full'])).max())
