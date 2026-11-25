import sys, json, time, numpy as np, math
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mylaw as M
P=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_predictions.json'))
FS={x['id']:x for x in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
c=P['curves'][0]; f=FS[c['id']]
cv=M.Curve(f['ainv'], f['N'])
FJ=cv.FJ(); FK=cv.FK()
res={}
for n in (10,20,40,80):
    tn,wn=M.nodes(n=n); Ln=M.Lt(f['N'],tn); wt=wn*Ln; wt/=wt.sum()
    res[n]=sum(w*(M.GUE+M.s2_from_F(FJ,FK,L)+cv.diag_corr(L)) for L,w in zip(Ln,wt))
for n in res: print(n, np.round(res[n]-res[80],6)[[3,7,11,13,15,17,19]])
# which part oscillates: the diag_corr part or the F part
for n in (10,80):
    tn,wn=M.nodes(n=n); Ln=M.Lt(f['N'],tn); wt=wn*Ln; wt/=wt.sum()
    a=sum(w*M.s2_from_F(FJ,FK,L) for L,w in zip(Ln,wt)); b=sum(w*cv.diag_corr(L) for L,w in zip(Ln,wt))
    print(n,'F part', np.round(a[[11,15,19]],6),'corr part', np.round(b[[11,15,19]],6))
