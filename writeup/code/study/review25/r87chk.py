import json, math, numpy as np, gmpy2
P = np.array([int(x) for x in open('/home/bas/praime/writeup/data/study/r17_primes.txt')], dtype=np.int64)
print(len(P), P[0]-10**18, P[-1]-10**18, 'lam', math.log(1e18+1e6), 480/math.log(1e18+1e6))
J = json.load(open('/home/bas/praime/writeup/data/study/r87.json'))
g = np.diff(P); rng = np.random.default_rng(0); QS=(3,5,7,11,13,17,19,23)
R = 400; acc = {(k,q):[] for k in range(1,7) for q in QS}
for rep in range(R):
    gs = rng.permutation(g); cs = np.concatenate([[0], np.cumsum(gs)])
    for k in range(1,7):
        d = cs[k:]-cs[:-k]; d = d[d<=480]
        for q in QS: acc[(k,q)].append(np.mean(d%q==0))
zl=[]; print('k q  law  data  indep  (law-indep)/se  (data-indep)/se')
sep = {}
for k in range(1,7):
    for q in QS:
        v = J['data'][str(k)][str(q)]; ind = np.mean(acc[(k,q)]); se = v['se']
        sep[(k,q)] = (J['theory'][str(k)][str(q)]-ind)/se
        if k>=2: print(k,q,'%.4f %.4f %.4f %+.2f %+.2f'%(J['theory'][str(k)][str(q)], v['p'], ind, sep[(k,q)], (v['p']-ind)/se))
for k in range(2,7):
    s=[sep[(k,q)] for q in QS]; print('k',k,'max |law-indep|/se %.2f'%max(map(abs,s)), 'sum sq %.1f'%sum(x*x for x in s))
print('total k>=2 sum (law-indep)^2/se^2 over 40: %.1f' % sum(sep[(k,q)]**2 for k in range(2,7) for q in QS))
# data vs indep chi2-like
