# True holomorphic newforms (trivial character, squarefree level, rational coefficients; coefficients from PARI) in the
# truncated forms: smallest eigenvalue of each form at the true level N (conductor N^k for Sym^k).  Usage: truth_hol.py <X>
import sys, numpy as np, cypari2
from maass import Problem
pari = cypari2.Pari(); X = float(sys.argv[1])
for N, k in [(1, 12), (1, 16), (11, 2), (5, 4), (3, 6), (2, 8), (2, 10), (14, 2), (15, 2), (7, 4)]:
    mf = pari.mfinit([N, k, 1], 0); B = pari.mfeigenbasis(mf)
    for F in B:
        co = pari.mfcoefs(F, 40)
        try: a = [int(co[n]) for n in range(41)]
        except Exception: continue
        P = Problem(float(k), 0, X, hol=True, J=3); two_a = 2*P.B.a; out = []
        for A, (name, _, kk, cf, _q) in zip(P.A, P.T):
            coef = []
            for p, m in P.pp:
                ap = a[p]/p**((k - 1)/2)
                if N % p == 0: w = np.sign(ap)**(kk*m)*p**(-kk*m/2)            # Steinberg: a_p = +-p^(-1/2)
                else: w = cf(m*np.arccos(ap/2))
                coef.append(w*np.log(p)*p**(-m/2))
            H = A + kk*np.log(N)*two_a*np.eye(P.B.n) - np.einsum('i,ijk->jk', np.array(coef), P.Spm)
            out.append('%s %+.2e' % (name, np.linalg.eigvalsh(H)[0]/(two_a*kk)))
        print('level %d weight %d X %g: %s   a_p/p^((k-1)/2): %s' % (N, k, X, '  '.join(out), ' '.join('%d:%+.4f' % (p, a[p]/p**((k - 1)/2)) for p in P.primes)), flush=True)
