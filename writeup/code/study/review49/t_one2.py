import sys, json, time, numpy as np, math
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mylaw as M
P=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_predictions.json'))
FS={x['id']:x for x in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
c=P['curves'][0]; f=FS[c['id']]
t0=time.time(); cv=M.Curve(f['ainv'], f['N']); print('init', time.time()-t0)
A2=cv._aot(pmax=20000); print('Aot 2e4 vs 1e5 max rel', np.abs(A2/cv.Aot-1).max())
FJ=cv.FJ(); FK=cv.FK()
for n in (10,12,20):
    tn,wn=M.nodes(n=n); Ln=M.Lt(f['N'],tn); wt=wn*Ln; wt/=wt.sum()
    law=sum(w*(M.GUE+M.s2_from_F(FJ,FK,L)+cv.diag_corr(L)) for L,w in zip(Ln,wt))
    print(n,'nodes: own-frozen', np.round(law-np.array(c['full']),6)[[0,3,7,11,15,19]])
# step h check
tn,wn=M.nodes(n=12); Ln=M.Lt(f['N'],tn); wt=wn*Ln; wt/=wt.sum()
lawh=sum(w*(M.GUE+M.s2_from_F(FJ,FK,L,h=0.001)+cv.diag_corr(L)) for L,w in zip(Ln,wt))
law=sum(w*(M.GUE+M.s2_from_F(FJ,FK,L)+cv.diag_corr(L)) for L,w in zip(Ln,wt))
print('h 0.001 vs 0.0025', np.abs(lawh-law).max())
