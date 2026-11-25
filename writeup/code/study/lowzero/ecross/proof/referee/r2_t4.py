"""Referee check 2 (own code): Theorem 4 (Voronoi identity) three ways, Gaussian-edged box, twist r.

  V_n   = sum_n lambda(n) n^{-1/2} int H(t) e^{i(2theta(t) - t ln(n r))} dt        (all integers n in the relevant range;
          I_n by the trapezoid rule, z^j recursion over the uniform grid)
  V_c   = 2pi * (1/2pi) int h(5/4+it) X^{-1}(5/4+it) L(5/4+it) r^{-3/4-it} dt          (PARI lfun on the 5/4 line)
  V_d   = w sum_m lambda(m) m^{-1/2} Hhat(ln(r/m)),  Hhat(xi) = int H e^{-i xi t} dt = What(xi) e^{-D^2 xi^2/4}
Theorem 4: V_n = V_d; contour form: V_c = V_n (per-n shift) and V_c = V_d (functional equation + shift to -1/4).
usage: r2_t4.py CID T1 T2 Delta r [margin]
"""
import os, sys, json, math, time
import numpy as np, mpmath as mp, cypari2
from scipy.special import erfc

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 7)))
cid = sys.argv[1]; T1, T2, D, r = map(float, sys.argv[2:6])
margin = float(sys.argv[6]) if len(sys.argv) > 6 else 4.0
C = {c['id']: c for c in json.load(open(REPO + '/writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}[cid]
N = C['N']
mp.mp.dps = 25
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6); pari.set_real_precision(38)
pari('E = ellinit(%s);' % str(C['ainv']))
w = int(pari('ellrootno(E)'))                     # root number from PARI, not from the json
PI = math.pi; t_start = time.time()

def theta(t):
    return float(t / 2 * mp.log(N) - t * mp.log(2 * mp.pi) + mp.im(mp.loggamma(1 + 1j * mp.mpf(t))))
thetav = np.vectorize(theta)
def erfdiff(a, b):
    a = np.asarray(a, complex); b = np.asarray(b, complex)
    return np.where((a.real + b.real) > 0, erfc(b) - erfc(a), erfc(-a) - erfc(-b))
def H(z):
    z = np.asarray(z, complex)
    return 0.5 * erfdiff((z - T1) / D, (z - T2) / D)
def Xinv(s):
    s = mp.mpc(s)
    return complex(mp.power(N, s - 0.5) * mp.power(2 * mp.pi, 1 - 2 * s) * mp.gamma(s + 0.5) / mp.gamma(1.5 - s))
lo, hi = T1 - 9 * D, T2 + 9 * D

# ---- V_d (dual side), m up to where e^{-D^2 ln^2(m/r)/4} < 1e-30
mmax = int(r * math.exp(2 * math.sqrt(70) / D)) + 2
am = np.array(pari('Vecsmall(ellan(E, %d))' % mmax), float); m = np.arange(1, mmax + 1, dtype=float)
xi = np.log(r / m)
What = np.where(np.abs(xi) < 1e-15, T2 - T1, (np.exp(-1j * xi * T2) - np.exp(-1j * xi * T1)) / (-1j * np.where(np.abs(xi) < 1e-15, 1, xi)))
Vd_terms = w * am / m * What * np.exp(-D * D * xi * xi / 4)          # lambda(m) m^{-1/2} = a_m / m
Vd = Vd_terms.sum()
# check one Hhat by quadrature (guards the closed form)
tq = np.arange(lo, hi, 1 / 64.0)
Hq = np.sum(H(tq) * np.exp(-1j * xi[1] * tq)) / 64.0
print('%s N=%d w=%+d (PARI) box [%g,%g] Delta=%g r=%g: dual side V_d = %.12f %+.12fi (m <= %d; largest dual terms: %s); Hhat(ln(r/2)) closed %.10f%+.10fi vs quadrature %.10f%+.10fi' % (
    cid, N, w, T1, T2, D, r, Vd.real, Vd.imag, mmax,
    ', '.join('m=%d:%.4f' % (k + 1, abs(Vd_terms[k])) for k in np.argsort(-np.abs(Vd_terms))[:4]),
    (What[1] * math.exp(-D * D * xi[1] ** 2 / 4)).real, (What[1] * math.exp(-D * D * xi[1] ** 2 / 4)).imag, Hq.real, Hq.imag))

# ---- V_c: contour integral on the 5/4 line with PARI L
pari('Lf = lfuncreate(E); Li = lfuninit(Lf, [1.75, 0.1, %d]);' % int(hi + 5))
step = 0.05
tt = np.arange(lo, hi + step / 2, step); vals = np.empty(len(tt), complex)
for j, t in enumerate(tt):
    Lv = complex(pari('lfun(Li, %r + I*%r)' % (1.75, float(t))))
    vals[j] = H(t - 0.75j) * Xinv(mp.mpc(1.25, t)) * Lv * r ** (-0.75) * np.exp(-1j * t * math.log(r))
wt = np.full(len(tt), step); wt[0] = wt[-1] = step / 2
Vc = np.sum(vals * wt)                                                    # = 2pi * (1/2pi i) int ds
print('V_c (5/4 line, PARI L) = %.12f %+.12fi   (%.0fs)' % (Vc.real, Vc.imag, time.time() - t_start))

# ---- V_n: sum over all integers n with ln(n r) within margin of [2theta'(T1), 2theta'(T2)]
def tthp(t):
    return float(mp.log(N) - 2 * mp.log(2 * mp.pi) + 2 * mp.re(mp.digamma(1 + 1j * mp.mpf(t))))
xl = math.exp(tthp(T1) - margin) / r; xh = math.exp(tthp(T2) + margin) / r
nmin = max(1, int(xl)); nmax = int(xh) + 1
a = np.array(pari('Vecsmall(ellan(E, %d))' % nmax), float)[nmin - 1:]
n = np.arange(nmin, nmax + 1, dtype=float)
coef = a / n                                                              # lambda(n) n^{-1/2}
del a
ht = 1 / 16.0
tg = np.arange(lo, hi + ht / 2, ht)
F = H(tg) * np.exp(2j * thetav(tg)) * ht; F[0] *= 0.5; F[-1] *= 0.5
Vn = 0j
for s in range(0, len(n), 400000):
    lnr = np.log(n[s:s + 400000] * r)
    z = np.exp(-1j * ht * lnr)
    acc = np.zeros(len(lnr), complex)
    for Fj in F[::-1]:
        acc = acc * z + Fj
    Vn += np.sum(coef[s:s + 400000] * acc * np.exp(-1j * tg[0] * lnr))
print('V_n (all integers %d..%d, %d terms) = %.12f %+.12fi   (%.0fs)' % (nmin, nmax, len(n), Vn.real, Vn.imag, time.time() - t_start))
print('|V_n - V_d| = %.2e   |V_c - V_n| = %.2e   |V_c - V_d| = %.2e   |V_n - (sign-flipped dual)| = %.2e' % (
    abs(Vn - Vd), abs(Vc - Vn), abs(Vc - Vd), abs(Vn + Vd)))
ok = abs(Vn - Vd) < 1e-8 and abs(Vc - Vd) < 1e-8
print('r2', 'PASS' if ok else 'FAIL', '(%.0fs)' % (time.time() - t_start))
