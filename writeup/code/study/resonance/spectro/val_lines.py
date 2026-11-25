"""val_lines.py LABEL: (1) bad-prime factors A_l(k) of the line r = l^k, measured exactly on integers sieved at the bad
prime l alone (weight l/(l-1) [l !| n]): A_l(k) = sum_n lambda(n) f_l(n) n^{-1/4} W(sqrt n) e(-2 sqrt(l^k n / N)) / MT(l^k);
(2) check of the line formula (lines_lib) on integers sieved at P = 10 (prime weight) for denominators with squares.
Integers n in [1e5, 1e7) (not primes).  Output badfac_<label>.json and a printed table."""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import numpy as np
import spectro_lib as S
import lines_lib as LL

label = sys.argv[1]; X0, X1 = 100000, 10_000_000
I = S.label_info(label); N = I['N']; v0, v1 = math.sqrt(X0), math.sqrt(X1)
bads = {'EC:11a1': (11, range(-2, 3)), 'EC:37a1': (37, range(-1, 2)), 'C23': (23, range(-1, 2)), 'A4': (2, range(-3, 5))}[label]
if label.startswith('EC:'):
    ptab, atab = S.ap_table(label[3:], 30_000_000); lamseg = lambda L, U: S.lam_ec_seg(L, U, ptab, atab, N)
elif label == 'A4': lamseg = lambda L, U: S.lam_gauss_seg(L, U, 4)
else: lamseg = lambda L, U: S.lam_c23_seg(L, U)
lb, ks = bads
RS2 = [Fr(1, 4), Fr(3, 4), Fr(9, 4), Fr(1, 9), Fr(2, 9), Fr(4, 9), Fr(1, 36), Fr(5, 4), Fr(1, 12), Fr(1, 8), Fr(1, 27), Fr(7, 9), Fr(1, 18)]
acc_bad = {k: 0j for k in ks}; acc2 = {r: 0j for r in RS2}; nz2 = {r: 0.0 for r in RS2}
t0 = time.time()
for L in range(X0, X1, 1_000_000):
    U = min(X1, L + 1_000_000)
    n = np.arange(L, U, dtype=np.int64); nf = n.astype(np.float64)
    w0 = lamseg(L, U) * nf ** -0.25 * S.hann(np.sqrt(nf), v0, v1)
    fb = (lb / (lb - 1.0)) * (n % lb != 0)
    for k in ks:
        r = float(lb) ** k
        acc_bad[k] += complex(np.dot(w0 * fb, np.exp(-1j * S.TWO_PI * 2 * np.sqrt(r * nf / N))))
    F = np.ones(U - L)
    for l in (2, 3, 5, 7): F *= (l / (l - 1.0)) * (n % l != 0)
    for r in RS2:
        acc2[r] += complex(np.dot(w0 * F, np.exp(-1j * S.TWO_PI * 2 * np.sqrt(float(r) * nf / N))))
        nz2[r] += float(np.dot(w0 * F, w0 * F))
bad = {}
for k in ks:
    MT = S.main_term_unit(label, v0, v1, float(lb) ** k); z = acc_bad[k] / MT
    bad['%d|%d' % (lb, k)] = z.real
    print('bad prime %d, r = %d^%d: A = %+.5f %+.5fi' % (lb, lb, k, z.real, z.imag))
json.dump(bad, open(LL.bad_factor_file(label), 'w'), indent=1)
# truncated prediction at P = 10: prod_{l <= 7} L_l(r) * prod_{l > 7, l | a} lambda(l^e) (b must be 7-smooth)
print('integers sieved at P = 10 (prime weight), lines with square denominators:')
for r in RS2:
    a, b = r.numerator, r.denominator
    pred = 1.0
    for l in (2, 3, 5, 7):
        lp, chi = LL.loc(label)[l]; ea, eb = LL.vval(a, l), LL.vval(b, l)
        if chi != 0.0: pred *= LL.L_good(label, l, ea, eb)
        else:
            key = '%d|%d' % (l, ea - eb); pred *= bad.get(key, float('nan'))
    for l, e in LL.factor(a).items():
        if l > 7: pred *= LL.lam_pk(label, l, e)
    MT = S.main_term_unit(label, v0, v1, float(r)); z = acc2[r] / MT; err = math.sqrt(nz2[r] / 2) / abs(MT)
    print('  r = %-5s read %+.4f %+.4fi   pred %+.4f   (random-model err %.4f)' % (r, z.real, z.imag, pred, err))
print('done %.0fs' % (time.time() - t0))
