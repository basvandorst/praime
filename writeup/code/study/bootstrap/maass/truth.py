# The true Maass forms in the truncated Weil forms: with the Hecke eigenvalues of the answer key (point masses at the
# true Satake angles), the smallest eigenvalue of each of the four forms, per unit conductor exponent.  Under GRH for
# L(Sym^k f) every number must be >= 0 (up to float64 and the basis).  Usage: truth.py <X> [rmax]
import sys, json, os, numpy as np
from maass import Problem
key = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass', 'lmfdb_key.json')))
X = float(sys.argv[1]); rmax = 20.0
for f in key:
    if f['r'] > rmax or not f['ap']: continue
    eps = 1 if f['parity'] == 'odd' else 0
    P = Problem(f['r'], eps, X, J=3, center=(eps == 1 and 'center' in sys.argv))
    th = {int(p): np.arccos(a/2) for p, a in f['ap'].items()}
    two_a = 2*P.B.a; out = []
    for A, (name, _, k, cf, _q) in zip(P.A, P.T):
        coef = np.array([cf(m*th[p])*np.log(p)*p**(-m/2) for p, m in P.pp])
        ev = np.linalg.eigvalsh(A - np.einsum('i,ijk->jk', coef, P.Spm))
        out.append('%s %+.2e' % (name, ev[0]/(two_a*k)))
    print('r %.8f %-4s X %g N %d  smallest eigenvalue: %s   max|a_p| %.3f' % (f['r'], f['parity'], X, P.N, '  '.join(out), max(abs(a) for a in f['ap'].values())), flush=True)
