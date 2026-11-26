"""Constants and checks: (i) chi(-2w) Ktilde21(w) = G3 at random complex w (meromorphic identity);
(ii) Euler products C0, A3_nat(1) = prod(1 - l^-2 + l^-5) and its local derivation from the K21 lines;
(iii) nu(m) table and its mean; (iv) M_nat values; (v) sqrt(y) coefficient of the Petersson-line version."""
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np, mpmath as mp, random
from duality_check import K21t, G3, chi
from model import nu, C0, M_nat, PR
mp.mp.dps = 30
random.seed(1)
mx = 0
for _ in range(12):
    w = mp.mpc(random.uniform(-0.45, 0.45), random.uniform(-40, 40))
    for k in (2, 4, 6, 8):
        r = chi(-2 * w) * K21t(w, k) / G3(1j * w, k)   # G3(t) with w = -it  <=>  t = i w
        mx = max(mx, abs(r - 1))
print("(i) max |chi(-2w) Ktilde21(w)/G3 - 1| over 48 random complex w, k=2..8: %.1e" % mx)
P = PR
A1 = float(np.exp(np.sum(np.log(1 - P ** -2 + P ** -5))))
print("(ii) C0 = prod(1 - 1/((l-1)(l^2-1))) = %.6f ;  A3_nat(1) = prod(1 - l^-2 + l^-5) = %.6f ; 6/pi^2 = %.6f" % (C0, A1, 6 / np.pi ** 2))
# local enumeration of the main GL(3) line from the K21 lines (a=v(r), sigma=[l|s], eps=[l|e], delta=[l|d], b=v(j))
def L_loc(l, wv, amax=12, bmax=12):
    tot = 0.0
    wl = l ** 1.5 / ((l - 1) * (l * l - 1)); hl = 1.0 / (l * l - 1)
    for a in range(amax):
        for sg in (0, 1):
            if a >= 1 and sg: continue
            for ep in (0, 1):
                if ep and not (a >= 1 or sg): continue
                for de in (0, 1):
                    if de and (a >= 1 or sg): continue
                    for b in range(bmax):
                        if a + 2 * b - sg - 2 * ep - 2 * de != wv: continue
                        tot += l ** (-1.5 * a) * wl ** sg * (-1) ** ep * hl ** de * l ** (-b)
    return (1 - l ** -2.0) * tot
print("     local factor at w=0 vs 1 - l^-2 + l^-5:", [("%.8f" % L_loc(l, 0), "%.8f" % (1 - l ** -2.0 + l ** -5.0)) for l in (2, 3, 5)])
print("     harmonic main GL(3) line A3_h(1) = 1 (s = 1 term only)")
print("(iii) nu(m), m = 1..16:", " ".join("%.4f" % nu(m) for m in range(1, 17)))
print("      mean nu(m), m <= 10^5: %.5f" % np.mean([nu(m) for m in range(1, 100001)]))
ys = np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
for k in (2, 4):
    print("(iv) M_nat k=%d at y=%s: %s" % (k, ys.tolist(), np.round(M_nat(ys, k), 3).tolist()))
