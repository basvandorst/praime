import numpy as np, json, sys
d = np.load('out/w_6000000000.npz')
J = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r135b_new1.json'))
for q,k in ((3,1),(3,2),(5,1),(7,1),(10,2)):
    for H in (10**4,10**5,10**6):
        w = d['win_%d_%d_%d'%(q,k,H)]; NB=len(w); mine = w.mean()*NB/(NB-1)
        theirs = J['win_%d_%d_%d'%(q,k,H)]['per_window'][0]
        print(q,k,H,'%.10f %.10f diff %.2e'%(mine,theirs,mine-theirs))
print(d['nprimes'], d['first_last'])
