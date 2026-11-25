"""Is there a common (curve-independent) t-structure in the genus-2 residual?  Window means (10 units) for independent halves
(round 151 vs round 152; random halves of the 55), their correlation, and the chi2 of window means vs their jackknife errors.
Also the same for Sym^3 (75) and the imaginary parts."""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
rows = load()
def windows(S, edges, part='real'):
    W = np.zeros((len(S), len(edges) - 1))
    for i, r in enumerate(S):
        g = r['g']; ph = np.exp(2j * theta(g, r['Q'], r['mu']))
        for j in range(len(edges) - 1):
            a, b = edges[j], edges[j + 1]
            m = (g > a) & (g <= b)
            v = ph[m].sum() + r['eps'] * r['c'] * (b - a) / (2 * math.pi)
            W[i, j] = v.real if part == 'real' else v.imag
    return W
for cls, edges in (('g2', np.arange(10, 151, 10)), ('sym3', np.arange(10, 51, 5))):
    S = [r for r in rows if r['cls'] == cls]
    for part in ('real', 'imag'):
        W = windows(S, edges, part)
        m = W.mean(0); se = W.std(0, ddof=1) / math.sqrt(len(S))
        chi2 = np.sum((m / se) ** 2)
        print('%s %s: window means %s' % (cls, part, np.round(m, 2)))
        print('      se %s ; chi2 vs 0 = %.1f on %d windows' % (np.round(se, 2), chi2, len(m)))
        # correlation between window means of independent halves (random splits)
        rng = np.random.default_rng(1)
        cors = []
        for k in range(2000):
            p = rng.permutation(len(S)); A = W[p[: len(S) // 2]].mean(0); B = W[p[len(S) // 2:]].mean(0)
            cors.append(np.corrcoef(A - A.mean(), B - B.mean())[0, 1])
        # null: shuffle window labels independently per curve (destroys common structure but keeps per-curve totals' distribution)
        cn = []
        for k in range(2000):
            Wp = np.array([w[rng.permutation(len(w))] for w in W]) if False else W * rng.choice([-1, 1], size=W.shape)
            p = rng.permutation(len(S)); A = Wp[p[: len(S) // 2]].mean(0); B = Wp[p[len(S) // 2:]].mean(0)
            cn.append(np.corrcoef(A - A.mean(), B - B.mean())[0, 1])
        print('      split-half correlation of window profiles: mean %.3f (null with random signs: mean %.3f, sd %.3f)' % (np.mean(cors), np.mean(cn), np.std(cn)))
        if cls == 'g2':
            r151 = np.array([r['rnd'] == 151 for r in S])
            A = W[r151].mean(0); B = W[~r151].mean(0)
            print('      round151 profile %s\n      round152 profile %s\n      corr %.3f' % (np.round(A, 2), np.round(B, 2), np.corrcoef(A, B)[0, 1]))
