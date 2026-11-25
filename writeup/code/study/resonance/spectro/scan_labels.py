"""scan_labels.py: labels a(p) for all primes p < 3e7 (DISCOVERY range), cached in labels_3e7.npz.
  l11 = a_p(11a)/sqrt p, l37 = a_p(37a)/sqrt p, c23 = cubic x^3 - x - 1 (2, -1, 0; 1 at 23), a4 = 2 cos 4 theta_p (0 if p = 3 mod 4,
  -1 at p = 2), rnd = random completely multiplicative +-1 (seed 20261125; on primes: iid +-1),
  prod = l11 * l37 (Rankin-Selberg 11a x 37a), sym2 = l11^2 - 1 (Sym^2 11a; 1/11 at p = 11)."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spectro_lib as S

X = 30_000_000
p11, a11 = S.ap_table('11a1', X)
p37, a37 = S.ap_table('37a1', X)
assert np.array_equal(p11, p37)
p = p11; pf = p.astype(np.float64)
l11 = a11 / np.sqrt(pf); l37 = a37 / np.sqrt(pf)
c23 = np.zeros(len(p)); a4 = np.zeros(len(p))
SEG = 2_000_000
for L in range(0, X, SEG):
    U = min(X, L + SEG); m = (p >= L) & (p < U); idx = p[m] - L
    c23[m] = S.lam_c23_seg(L, U)[idx]
    a4[m] = S.lam_gauss_seg(L, U, 4)[idx]
rng = np.random.default_rng(20261125)
rnd = rng.choice([-1.0, 1.0], size=len(p))
sym2 = l11 * l11 - 1.0; sym2[p == 11] = 1.0 / 11
np.savez(os.path.join(S.HERE, 'labels_3e7.npz'), p=p, l11=l11, l37=l37, c23=c23, a4=a4, rnd=rnd, prod=l11 * l37, sym2=sym2)
# sanity: a few values
for q in (2, 3, 5, 7, 11, 13, 23, 29, 37, 59):
    k = np.searchsorted(p, q)
    print(q, 'l11 %.4f l37 %.4f c23 %g a4 %.4f' % (l11[k], l37[k], c23[k], a4[k]))
print('second moments: l11 %.4f l37 %.4f c23 %.4f a4 %.4f rnd %.4f prod %.4f sym2 %.4f' % tuple(np.mean(x ** 2) for x in (l11, l37, c23, a4, rnd, l11 * l37, sym2)))
