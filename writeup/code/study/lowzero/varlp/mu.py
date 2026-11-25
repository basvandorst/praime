"""mu.py: the Palm mean mu(L) = int ln|r| T2(r) dr of the base (q = 1 arithmetic) at fixed L, which enters the
variance through -mu^2 (derivation.md sect. 4.4):
    mu = mu_GUE + (gamma + C1)/L + mu_off(L),   mu_GUE = gamma + ln 2pi - 1,
    C1 = sum_p sum_{m>=2} (1 - 1/m) ln p / p^m           (Mertens: the diagonal part, exact),
    mu_off = int ln|r| (2/L^2) Re[e^{2 pi i r} rego(eps)] dr  (numerical, damped; converges fast).
Check: the direct damped integral of the whole C(r) (slow, ~ ln E0/E0) approaches the same number.
usage: mu.py OUT.json
"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
from oned import rnodes
EG = vlib.EG; TWO_PI = vlib.TWO_PI
P = vlib.primes_upto(10 ** 7).astype(float)
C1 = float(sum((np.log(P) * (1 - 1 / m) * P ** (-float(m))).sum() for m in range(2, 70)))
print('C1 = %.8f' % C1)
out = {'C1': C1, 'mu_GUE': EG + math.log(TWO_PI) - 1, 'rows': []}
for L in (7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0):
    Rmax = 150.0 * L / TWO_PI; r, w = rnodes(Rmax); e = TWO_PI * r / L; lr = np.log(r)
    rego = vlib.X('Z2m1', e) * vlib.X('Phi', e) + vlib.X('Phim1', e)
    off = (2 / L ** 2) * (np.exp(1j * TWO_PI * r) * rego).real
    C = vlib.C_of_r(L, (), r)
    vals_off = [2 * np.sum(w * lr * off * np.exp(-(e / E0) ** 2)) for E0 in (15., 30., 45.)]
    vals_tot = [2 * np.sum(w * lr * C * np.exp(-(e / E0) ** 2)) for E0 in (15., 30., 45.)]
    mu_ar = (EG + C1) / L + vals_off[1]
    out['rows'].append(dict(L=L, diag=(EG + C1) / L, off=vals_off, mu_arith=mu_ar, direct_damped=vals_tot, mu=out['mu_GUE'] + mu_ar))
    print('L=%5.2f  diag (gamma+C1)/L %.5f  off (E0=15/30/45) %s  -> mu_arith %.5f, mu %.5f | direct damped C: %s' % (
        L, (EG + C1) / L, ' '.join('%+.6f' % x for x in vals_off), mu_ar, out['mu_GUE'] + mu_ar, ' '.join('%.5f' % x for x in vals_tot)))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
