"""Round 114 (preregister/round114.txt): S4 quartic fields |disc| <= 20000, Artin L = zeta_K/zeta: first zero, L(1), first
prime with a root, and the Hadamard sum S; Spearman correlations and the S-prediction of the first zero. usage: r111.py"""
import json, math, os
import numpy as np, cypari2
from scipy.stats import spearmanr
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
rows = []
for sig in (0, 1, 2):
    for f in pari('nflist("S4",[1,20000],%d)' % sig):
        disc = int(pari.nfdisc(f)); N = abs(disc)
        A = pari('lfundiv(lfuncreate(nfinit(%s)), lfuncreate(1))' % f)
        T = 6; z = []
        while not z: z = [float(t) for t in pari.lfunzeros(A, T) if float(t) > 0]; T *= 2
        L1 = float(pari.lfun(A, 1)); dL1 = float(pari.lfun(A, 1, 1))
        mu = [int(m) for m in pari.lfunparams(A)[2]]
        S = 2 * dL1 / L1 + math.log(N / math.pi ** len(mu)) + sum(float(pari.psi((1 + m) / 2.0)) for m in mu)
        P = next(int(p) for p in pari.primes(300) if N % int(p) and int(pari.polrootsmod(f, p).length()) > 0) if any(N % int(p) and int(pari.polrootsmod(f, p).length()) > 0 for p in pari.primes(300)) else 2000
        rows.append({'f': str(f), 'disc': disc, 'sig': sig, 'gamma1': min(z), 'L1': L1, 'S': S, 'P': P, 'g': min(z) * math.log(N) / (2 * math.pi)})
out = {'n': len(rows), 'tests': {}}
for sig in (0, 1, 2):
    R = [r for r in rows if r['sig'] == sig]
    lnN = np.array([math.log(abs(r['disc'])) for r in R]); S = np.array([r['S'] for r in R]); g1 = np.array([r['gamma1'] for r in R])
    pred = np.empty(len(R))
    for i in range(len(R)):
        m = np.abs(lnN - lnN[i]) < 0.25
        lhs = S[i] - np.median(S[m]) + 2 / (0.25 + np.median(g1[m]) ** 2)
        pred[i] = math.sqrt(max(2 / lhs - 0.25, 0)) if lhs > 0 else np.nan
    ok = ~np.isnan(pred)
    t = {'n': len(R), 'spearman_lnL1_g': spearmanr([math.log(r['L1']) for r in R], [r['g'] for r in R])[0],
         'spearman_lnP_g': spearmanr([math.log(r['P']) for r in R], [r['g'] for r in R])[0],
         'spearman_pred_gamma1': spearmanr(pred[ok], g1[ok])[0], 'n_no_inversion': int((~ok).sum()),
         'pearson_S_inv': float(np.corrcoef(S, 2 / (0.25 + g1 ** 2))[0, 1]),
         'lowest5': sorted([(round(r['gamma1'], 4), r['disc'], r['P'], round(r['L1'], 3), r['f']) for r in R])[:5]}
    out['tests']['s=%d' % sig] = t
    print('signature', sig, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in t.items()})
out['rows'] = rows
json.dump(out, open(os.path.join(H, 'r114.json'), 'w'), indent=0)
