import numpy as np, json, math, glob
WIN = [6*10**9 + i*10**8 for i in range(24)]
D = [np.load('out/w_%d.npz' % x) for x in WIN]
HS = (10**4, 2*10**4, 5*10**4, 10**5, 2*10**5, 5*10**5, 10**6)
R = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r137.json'))['cells']
out = {}
for q, k in ((3,1),(5,1),(7,1),(3,2),(10,2),(10,1)):
    for H in HS:
        key = 'win_%d_%d_%d' % (q,k,H)
        pw = []; pe = []
        for d in D:
            w = d[key]; NB = len(w); f = NB/(NB-1)
            pw.append(w.mean()*f); pe.append(w.std()*f/math.sqrt(NB))
        pw = np.array(pw); pe = np.array(pe)
        m = pw.mean(); sb = math.sqrt((pe**2).sum())/24; ss = pw.std(ddof=1)/math.sqrt(24)
        out[key] = dict(mean=m, sb=sb, ss=ss, pw=pw.tolist())
        rk = '%d_%d_%d' % (q,k,H)
        if rk in R:
            r = R[rk]
            print('%-12s mine %.6f (sb %.5f ss %.5f) | r137 %.6f (sb %.5f ss %.5f)  diff %.1e' % (rk, m, sb, ss, r['meas'], r['sig_block'], r['sig_scatter'], m - r['meas']))
        elif H in (2*10**4, 5*10**4, 2*10**5, 5*10**5) or k == 1 and q == 10:
            print('%-12s mine %.6f (sb %.5f ss %.5f)' % (rk, m, sb, ss))
json.dump(out, open('agg.json','w'))
