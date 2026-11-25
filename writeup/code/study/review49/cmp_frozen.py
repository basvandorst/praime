import json, glob, numpy as np
P={c['id']:c for c in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_predictions.json'))['curves']}
D=[];K=[];names=[]
for f in sorted(glob.glob('lawout/*.json')):
    R=json.load(open(f)); c=P[R['id']]
    d=np.array(R['full']['full'])-np.array(c['full']); D.append(d)
    k=np.array(R['full']['deg1'])-np.array(c['deg1']); K.append(k)
    names.append(R['id'])
D=np.array(D); K=np.array(K)
print(len(D),'curves; own full - frozen full: max|.| per l:'); print(np.round(np.abs(D).max(0),5))
print('mean', np.round(D.mean(0),5))
print('own deg1 - frozen deg1 max', np.round(np.abs(K).max(0),4))
for nm,d in zip(names,D):
    if np.abs(d).max()>3e-4: print(nm, np.round(d,5))
