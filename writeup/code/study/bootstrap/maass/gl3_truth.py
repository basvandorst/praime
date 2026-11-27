# True SL(3,Z) forms (answer key gl3_key.json) in the two forms: smallest eigenvalues, which must be >= 0 under GRH.
# Also a check of the Hermitian engine against engine.py on a self-dual datum.  Usage: gl3_truth.py <X> [nforms]
import sys, json, os, numpy as np
from gl3 import GL3, HBasis
from engine import Basis
X = float(sys.argv[1]); nf = int(sys.argv[2]) if len(sys.argv) > 2 else 8
B1 = Basis(X, 16); B2 = HBasis(X, 16, 0.0)
print('engine check (zeta Gamma factor, no pole): max difference of the two archimedean matrices %.1e' % np.abs(B1.arch([0]) - B2.arch([0.0])).max())
print('engine check (Gamma_R(s + 3i) Gamma_R(s - 3i)): %.1e' % np.abs(B1.arch([3j, -3j]) - B2.arch([3.0, -3.0])).max())
key = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass', 'gl3_key.json')))
forms = sorted(key['forms'], key=lambda f: f['laplace'])[:nf]
for f in forms:
    l1, l2 = f['lambda']; P = GL3(l1, l2, X, ng=6); ang = []
    for p in P.primes:
        ap = complex(*f['A_p1'][str(p)])
        al = np.roots([1, -ap, np.conj(ap), -1])                              # Satake parameters: roots of x^3 - a x^2 + conj(a) x - 1
        ang.append(np.angle(al))
    dev = max(abs(abs(np.exp(1j*np.array(ang))) - 1).max(), 0)
    v = P.point(ang)
    vd = P.point([-np.array(x) for x in ang])                                 # the dual's coefficients with the same Gamma factor (must fail if orientation matters)
    print('(%.5f, %.5f) Laplace %.2f X %g N %d: F %+.2e  Ad %+.2e   | with conjugated a_p: F %+.2e   max||alpha|-1| %.1e' % (
        l1, l2, f['laplace'], X, P.N, v[0], v[1], vd[0], max(abs(abs(np.roots([1, -complex(*f['A_p1'][str(p)]), np.conj(complex(*f['A_p1'][str(p)])), -1])) - 1).max() for p in P.primes)), flush=True)
