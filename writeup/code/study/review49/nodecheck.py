import sys, json, glob, numpy as np
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
P={c['id']:c for c in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_predictions.json'))['curves']}
rows=[]
for f in glob.glob('lawout/*.json'):
    r=json.load(open(f)); d=np.array(r['full']['full'])-np.array(P[r['id']]['full']); rows.append((abs(d[19]),r['id'],r['N'],d[[15,17,19]]))
rows.sort(reverse=True)
for x in rows[:6]: print(x[1],x[2],np.round(x[3],4))
import mylaw as M
FS={x['id']:x for x in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
f=FS[rows[0][1]]; cv=M.Curve(f['ainv'],f['N']); FJ=cv.FJ(); FK=cv.FK()
res={}
for n in (10,24,48,96):
    tn,wn=M.nodes(n=n); Ln=M.Lt(f['N'],tn); wt=wn*Ln; wt/=wt.sum()
    res[n]=sum(w*(M.GUE+M.s2_from_F(FJ,FK,L)+cv.diag_corr(L)) for L,w in zip(Ln,wt))
for n in res: print(n,'nodes minus 96:',np.round(res[n]-res[96],5)[[11,15,17,19]], ' minus frozen:', np.round(res[n]-np.array(P[f['id']]['full']),5)[[15,17,19]])
