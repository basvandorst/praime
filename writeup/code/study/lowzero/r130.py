"""Round 130 (preregister/round130.txt): F5 quintic fields |disc| <= 1e7 (list in f5list.txt from PARI nflist("F5",[1,10^7])),
degree-4 Artin L = zeta_K/zeta: first zero predicted from S. usage: r130.py"""
import os, json, math
import numpy as np, cypari2
from scipy.stats import spearmanr
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in [l.strip() for l in open(os.path.join(H, 'f5list.txt')) if l.strip()]:
    A = pari('lfundiv(lfuncreate(nfinit(%s)), lfuncreate(1))' % f)
    N = int(pari.lfunparams(A)[0]); mu = [int(m) for m in pari.lfunparams(A)[2]]
    T = 2.0; z = []
    while not z: z = [float(t) for t in pari.lfunzeros(A, T) if float(t) > 0]; T *= 2
    L1 = float(pari.lfun(A, 1)); dL1 = float(pari.lfun(A, 1, 1))
    S = 2 * dL1 / L1 + math.log(N / math.pi ** len(mu)) + sum(float(pari.psi((1 + m) / 2.0)) for m in mu)
    rows.append({'f': f, 'N': N, 'mu': mu, 'gamma1': min(z), 'L1': L1, 'S': S, 'g': min(z) * math.log(N) / (2 * math.pi)})
lnN = np.array([math.log(r['N']) for r in rows]); S = np.array([r['S'] for r in rows]); g1 = np.array([r['gamma1'] for r in rows])
pred = np.empty(len(rows))
for i in range(len(rows)):
    m = np.abs(lnN - lnN[i]) < 0.5
    lhs = S[i] - np.median(S[m]) + 2 / (0.25 + np.median(g1[m]) ** 2)
    pred[i] = math.sqrt(max(2 / lhs - 0.25, 0)) if lhs > 0 else np.nan
ok = ~np.isnan(pred)
r1 = spearmanr(pred[ok], g1[ok])[0]; r2 = spearmanr([math.log(r['L1']) for r in rows], [r['g'] for r in rows])[0]
print('n %d (no inversion %d); Spearman pred-true %.3f; Spearman lnL1-g %.3f; readings %s %s' % (len(rows), int((~ok).sum()), r1, r2, r1 >= 0.9, r2 >= 0.4))
json.dump({'rows': rows, 'spearman_pred': float(r1), 'spearman_lnL1_g': float(r2), 'readings': [bool(r1 >= 0.9), bool(r2 >= 0.4)]}, open(os.path.join(H, 'r130.json'), 'w'), indent=0)
