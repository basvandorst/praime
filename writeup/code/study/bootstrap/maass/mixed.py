# Exploration (not registered): the odd Maass problem with the central zero, plus Rankin-Selberg constraints against
# forms already found.  For a known form g (parameter r_g, parity e_g, Hecke eigenvalues from the key) the L-function
# L(f x g) has Gamma factor Gamma_R(s + d +- i(r + r_g)) Gamma_R(s + d +- i(r - r_g)), d = (eps + e_g) mod 2, conductor 1,
# no pole when f != g, and weight (2 cos m theta_p)(2 cos m theta_p(g)) log p at p^m: linear in the unknown measure of f.
# Usage: mixed.py <X> <r,r,...> <number of known forms used: 0..4>
import sys, json, os, time, numpy as np
from maass import Problem, types
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
key = [f for f in json.load(open(os.path.join(D, 'lmfdb_key.json'))) if f['ap']]
X = float(sys.argv[1]); nk = int(sys.argv[3]); eps = 1
for r in [float(v) for v in sys.argv[2].split(',')]:
    T = types(r, eps, 2)
    for g in key[:nk]:
        if abs(g['r'] - r) < 1e-6: continue                                  # f x f has a pole: not used
        e = 1 if g['parity'] == 'odd' else 0; d = (eps + e) % 2; rg = g['r']
        thg = {int(p): np.arccos(a/2) for p, a in g['ap'].items()}
        fn = lambda ph: ph
        fn.pmp = (lambda thg: lambda p, m, th: 2*np.cos(m*th)*2*np.cos(m*thg[p]))(thg)
        T.append(('x%.2f' % rg, [d + 1j*(r + rg), d - 1j*(r + rg), d + 1j*(r - rg), d - 1j*(r - rg)], 2, fn))
    t0 = time.time(); P = Problem(r, eps, X, T=T, center=True, N=max(16, int(np.ceil(np.log(X)/2*(2*r + 12 + (key[nk - 1]['r'] if nk else 0))/np.pi))))
    lam, t, rd = P.solve()
    act = [x[0] for x, v in zip(P.T, P.mineigs(P.w)) if v < t + 1e-4]
    print('X %g r %.4f, %d known forms: t* = %+.5f (%d rounds, %.0fs) active %s; a_p %s' % (X, r, nk, t, rd, time.time() - t0, act, P.summary()), flush=True)
