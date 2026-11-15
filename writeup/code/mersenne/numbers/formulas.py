"""Direct formulas a_n = g(n) and their description length.

1. Least squares of ln a_n on n (a_n ~ A*B^n), and of ln a_n ~ alpha + beta*n with an offset C fitted.
2. Exact-hit maximum of two-parameter families, on the real sequence and on every null sequence:
     F1 floor(c * r^n)       : ln a_n   <= u + n v      < ln(a_n + 1)
     F2 floor(theta^(r^n))   : lnln a_n <= u + n v      < lnln(a_n + 1)       (Mills type, u = lnln theta, v = ln r)
     F3 floor(c * n^k)       : ln a_n   <= u + ln(n) v  < ln(a_n + 1)
     F4 round(c * r^n)       : ln(a_n - 1/2) <= u + n v < ln(a_n + 1/2)
   Each term is a strip in the (u, v) plane; the best point lies next to a vertex of the strip
   boundaries, so every pairwise boundary intersection is tried with four tiny nudges and the
   strips containing it are counted exactly (half-open test). Any two terms can always be hit.
3. Description length: raw bits sum(log2 a_n); Newton forward differences (integers, exact
   interpolating polynomial of degree 51); the Mills/Champernowne-style constant.
"""
import numpy as np, math, sys
from fractions import Fraction
D = np.load('seqs.npz')
real = D['real'].astype(np.int64)
L = len(real)
n = np.arange(1, L + 1, dtype=np.float64)
out = []
def pr(*a):
    s = ' '.join(str(x) for x in a); out.append(s); print(s)

# 1. least squares
y = np.log(real.astype(np.float64))
beta, alpha = np.polyfit(n, y, 1)
fit = np.exp(alpha + beta * n)
pr(f'LSQ ln a_n = {alpha:.4f} + {beta:.5f} n  ->  a_n ~ {math.exp(alpha):.4f} * {math.exp(beta):.5f}^n')
pr(f'   Wagstaff ratio 2^(1/e^gamma) = {2 ** math.exp(-0.5772156649015329):.5f}; empirical rate lambda = 51/ln(a52/a1) -> ratio e^(1/lambda) = {math.exp(math.log(real[-1] / 2) / 51):.5f}')
r = np.rint(fit).astype(np.int64)
pr('   exact hits round(A B^n):', int((r == real).sum()), ' within 1%:', int((np.abs(fit - real) <= 0.01 * real).sum()),
   ' median |ln error| = %.3f, max = %.3f' % (np.median(np.abs(y - np.log(fit))), np.max(np.abs(y - np.log(fit)))))
# with offset C: grid on C, least squares on ln(a_n - C)
bestC = None
for C in np.linspace(-1.9, 1.9, 39):
    yy = np.log(real - C); b2, a2 = np.polyfit(n, yy, 1); f2 = np.exp(a2 + b2 * n) + C
    h = int((np.rint(f2) == real).sum())
    if bestC is None or h > bestC[0]:
        bestC = (h, C, a2, b2)
pr(f'   round(A B^n + C), C on a grid of 39 values: best exact hits {bestC[0]} (C = {bestC[1]:.2f})')


# 2. strip arrangements
def strips(a, fam):
    a = a.astype(np.float64)
    if fam == 'F1': lo, hi, phi = np.log(a), np.log(a + 1), np.arange(1, len(a) + 1.0)
    if fam == 'F2': lo, hi, phi = np.log(np.log(a)), np.log(np.log(a + 1)), np.arange(1, len(a) + 1.0)
    if fam == 'F3': lo, hi, phi = np.log(a), np.log(a + 1), np.log(np.arange(1, len(a) + 1.0))
    if fam == 'F4': lo, hi, phi = np.log(a - .5), np.log(a + .5), np.arange(1, len(a) + 1.0)
    return lo, hi, phi


def best_hits(a, fam):
    lo, hi, phi = strips(a, fam)
    m = len(a)
    I, J = np.triu_indices(m, 1)
    pts = []
    for bi in (lo, hi):
        for bj in (lo, hi):
            v = (bj[J] - bi[I]) / (phi[J] - phi[I])
            u = bi[I] - phi[I] * v
            pts.append((u, v))
    U = np.concatenate([p[0] for p in pts]); V = np.concatenate([p[1] for p in pts])
    best, arg = 0, None
    global OOS
    OOS = (-1, 0)
    for du in (-1, 1):
        for dv in (-1, 1):
            # nudges well below the narrowest strip (~1/(a ln a)) but above float resolution
            uu = U + du * 1e-13 * np.maximum(1, np.abs(U)); vv = V + dv * 1e-13 * np.maximum(1, np.abs(V)) / 60
            val = uu[:, None] + phi[None, :] * vv[:, None]
            inside = (val >= lo[None, :]) & (val < hi[None, :])
            c = inside.sum(1)
            ci, co = inside[:, :25].sum(1), inside[:, 25:].sum(1)
            mi = int(ci.max())
            mo = int(co[ci == mi].max())
            if mi > OOS[0] or (mi == OOS[0] and mo > OOS[1]):
                OOS = (mi, mo)
            k = int(c.argmax())
            if c[k] > best:
                best, arg = int(c[k]), (uu[k], vv[k])
    return best, arg


def exact_check(a, fam, uv):  # recompute hits with the formula itself
    u, v = uv
    nn = np.arange(1, len(a) + 1, dtype=np.float64)
    if fam == 'F1': g = np.floor(np.exp(u + nn * v))
    if fam == 'F2': g = np.floor(np.exp(np.exp(u + nn * v)))
    if fam == 'F3': g = np.floor(np.exp(u) * nn ** v)
    if fam == 'F4': g = np.rint(np.exp(u + nn * v))
    return [int(i + 1) for i in np.nonzero(g == a)[0]]


nulls = {'A': D['nullA'], 'B': D['nullB']}
for fam, desc in (('F1', 'floor(c*r^n)'), ('F2', 'floor(theta^(r^n))'), ('F3', 'floor(c*n^k)'), ('F4', 'round(c*r^n)')):
    b, uv = best_hits(real, fam)
    idx = exact_check(real, fam, uv)
    if fam in ('F1', 'F4'): par = f'c = {math.exp(uv[0]):.6g}, r = {math.exp(uv[1]):.8g}'
    if fam == 'F2': par = f'theta = {math.exp(math.exp(uv[0])):.8g}, r = {math.exp(uv[1]):.8g}'
    if fam == 'F3': par = f'c = {math.exp(uv[0]):.6g}, k = {uv[1]:.8g}'
    line = f'{desc:20s} real: best {b} exact hits ({par}); recomputed hits at n = {idx}'
    b_r, uv_r = best_hits(real, fam); oos_real = OOS
    line += f' | OOS (fit on a_1..a_25): best in-sample {oos_real[0]}, max hits among a_26..a_52 over the tied fits {oos_real[1]}'
    for nm, Mx in nulls.items():
        bs, oo = [], []
        for s_ in Mx:
            bs.append(best_hits(s_, fam)[0]); oo.append(OOS[1])
        bs = np.array(bs); oo = np.array(oo)
        line += f' | null{nm} OOS max-over-ties mean {oo.mean():.2f}, share >= real {(1 + (oo >= oos_real[1]).sum()) / (1 + len(oo)):.3f}'
        pv = (1 + (bs >= b).sum()) / (1 + len(bs))
        line += f' | null{nm}: median {np.median(bs):.0f}, 95% {np.quantile(bs, .95):.0f}, max {bs.max()}, p = {pv:.3f}'
    pr(line)

# 3. description length
raw = sum(math.log2(int(v)) for v in real)
pr(f'\nraw description: sum log2 a_n = {raw:.0f} bits (Elias-style + 2 log2 log2 per term ~ {raw + sum(2 * math.log2(max(2, math.log2(int(v)))) for v in real):.0f})')
d = [int(v) for v in real]; newton = []
while d:
    newton.append(d[0]); d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
nb = sum(math.log2(abs(c)) + 1 for c in newton if c)
pr(f'Newton interpolation a_n = sum_k D^k a_1 * C(n-1, k): 52 integer coefficients, {nb:.0f} bits (largest |D^k| ~ 2^{max(math.log2(abs(c)) for c in newton if c):.0f})')
# monomial coefficients of the interpolating polynomial, exact rationals
from sympy import interpolate, symbols, Poly, Rational
xs = symbols('x')
poly = Poly(interpolate([(i + 1, int(v)) for i, v in enumerate(real)], xs), xs)
cb = sum(math.log2(abs(c.p)) + math.log2(c.q) + 1 for c in poly.all_coeffs() if c != 0)
pr(f'monomial form of the same polynomial: degree {poly.degree()}, rational coefficients totalling {cb:.0f} bits')
pr('Mills/Champernowne constant theta = 0.2 3 5 7 13 ... (each a_n written with a length prefix) encodes the list in '
   f'>= {raw:.0f} bits: it describes, it does not compress')
# MDL of the exponential fit: code ln a_n residuals with N(0, s^2) to integer precision
res = y - (alpha + beta * n); s = res.std(ddof=2)
# P(a_n) ~ density of ln a at y times width 1/a_n
bits_fit = sum(-math.log2(math.exp(-ri * ri / (2 * s * s)) / (s * math.sqrt(2 * math.pi)) / int(v)) for ri, v in zip(res, real)) + 2 * 32
pr(f'exponential model + Gaussian residuals in ln a (s = {s:.3f}): {bits_fit:.0f} bits incl. 64 for the 2 constants, vs raw {raw:.0f}')
lam = 51 / math.log(real[-1] / 2)
bits_pois = 1 + sum(-math.log2(lam * (int(real[i]) / int(real[i + 1])) ** lam / int(real[i + 1])) for i in range(51))
pr(f'Poisson-in-ln model (rate {lam:.3f}, 1 constant): {bits_pois:.0f} bits')
open('formulas.txt', 'w').write('\n'.join(out) + '\n')
