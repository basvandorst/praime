import json, math, numpy as np
r = np.load('rlaw.npz'); out = r['out']; lams = r['lams']
u = np.linspace(lams[0], lams[-1], 401); wt = np.exp(u) / u ** 2
cnt = json.load(open('counts.json'))
T=np.load('twins.npy'); W=T[T<2*10**9]; n=len(W)
for k in (1,2,3):
    m = np.sum([np.interp(u, lams, out[k-1].sum(1))*wt])/wt.sum()
    N = n - k; c = sum(cnt[f'{k}_in'].values()); f = c/N; se = math.sqrt(f*(1-f)/N)
    print(k, 'data frac <=1200 %.5f law %.5f z(naive) %.2f' % (f, m, (f-m)/se))
# shuffled-gap comparator from the data's own gaps (independence with the true marginal)
g = np.diff(W); rng = np.random.default_rng(0)
BANDS = [(lo, lo + 24) for lo in range(6, 600, 30)] + [(lo, lo + 54) for lo in range(606, 1200, 60)]
for k in (2,3):
    real = W[k:]-W[:-k]
    sh = np.zeros(len(BANDS))
    for rep in range(5):
        gs = rng.permutation(g); cs = np.concatenate([[0], np.cumsum(gs)]); d = cs[k:]-cs[:-k]
        sh += np.array([np.sum((d>=lo)&(d<=hi)) for lo,hi in BANDS])/5
    ob = np.array([np.sum((real>=lo)&(real<=hi)) for lo,hi in BANDS])
    sh *= ob.sum()/sh.sum()
    print(k, 'shuffled-gap chi2 %.0f' % np.sum((ob-sh)**2/sh))
