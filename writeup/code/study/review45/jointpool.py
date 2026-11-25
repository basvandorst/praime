import numpy as np, json, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from jointfit import joint, ols, PR
P = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/triple/predictions_kappa3.json'))
def arr(T):
    R = json.load(open(f__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/r138_{T}.json'))
    ex = [e[0] for e in R['excluded']]; Q = P[T]['D']; keep = [i for i, D in enumerate(Q) if D not in ex]
    a = {k: np.array(P[T][k])[keep] for k in ('full', 'Q0', 'diagonly')}
    a['k3'] = np.array(R['kappa3']); a['X'] = np.array([[1.] + [float(abs(Q[i]) % p == 0) for p in PR] for i in keep])
    assert len(a['k3']) == len(keep)
    return a
rng = np.random.default_rng(3); res = {}
for T in ('Fn', 'R120', 'B', 'Dn', 'C'):
    a = arr(T); n = len(a['k3']); g, s, rho, chi, dof = joint(a)
    bs = []
    for _ in range(400):
        ii = rng.integers(0, n, n); bs.append(joint({k: v[ii] for k, v in a.items()})[0])
    bs = np.array(bs); bse = bs.std(0, ddof=1); res[T] = (g, bse)
    print('%-5s n=%d joint: gamma_Q %.3f (model se %.3f, boot se %.3f), gamma_off %.3f (model %.3f, boot %.3f), corr %.2f' % (T, n, g[0], s[0], bse[0], g[1], s[1], bse[1], rho))
for lab, ts in (('all five', ('Fn', 'R120', 'B', 'Dn', 'C')), ('test O', ('R120', 'B', 'Dn', 'C'))):
    for k, nm in ((0, 'gamma_Q'), (1, 'gamma_off')):
        w = np.array([1/res[t][1][k]**2 for t in ts]); m = sum(wi*res[t][0][k] for wi, t in zip(w, ts))/w.sum()
        print('pooled %-8s %-9s (joint, boot se) %.4f +- %.4f' % (lab, nm, m, 1/np.sqrt(w.sum())))
