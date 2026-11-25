# rh_check.py -- checks every numeric claim and exponent computation of rh_proof.md.
# Self-contained (numpy, sympy, mpmath). Run:
#   (ulimit -v 1200000; timeout 600 nice /home/bas/praime/.venv/bin/python -I -B rh_check.py)
# Sections:
#   A  exponent of the large-values analysis as a function of kappa (Section 6 of the note)
#   B  mean-value lemma C' (Lemma 5): exact left side vs the proved bound, incl. the counterexample
#   C  divisor-moment bound with the content weight (Lemma 6)
#   D  lattice points in polar boxes (Lemma 2)
#   E  angular cut-off functions (Lemma 1)
#   F  Gaussian norm weight: Mellin transform and tail bounds (Lemma 3)
#   G  Heath-Brown identity on ideals of Z[i], sharp and smoothed Moebius cut (Lemma 4)
#   H  characters mod q, units, infinity types (Section 2)
#   I  Gaussian primes in residue classes mod 2+i and mod 5 in sectors of width c x^(-1/4), x = 1e7
#   J  parameter bookkeeping: an admissible choice of (theta, eps, Delta, J) for kappa = 0.24
#   K  large sieve (Gaussian-majorant form) and Gallagher's inequality on random data
import math, cmath, itertools, random
import numpy as np
import sympy as sp
import mpmath as mp

NPASS = 0
NFAIL = 0
def check(name, ok, info=""):
    global NPASS, NFAIL
    if ok:
        NPASS += 1
    else:
        NFAIL += 1
    print(("PASS " if ok else "FAIL ") + name + (("  [" + info + "]") if info else ""), flush=True)

E = math.e
SQPI = math.sqrt(math.pi)

# ---------------------------------------------------------------------------------------------
print("== A. exponent of the large-values analysis")
# Notation: M = family size (|m| ~ M), N = M^beta = length of the factor, sigma = large-value level.
# Available bounds for R (number of large values), as exponents of M (logs dropped):
#   trivial            R << M                                        -> 1
#   C' with power g    R << (M + N^g) N^{g(1-2 sigma)}               -> max(1 + g b (1-2s), 2 g b (1-s))
#   pointwise (long)   R = 0 if sigma > 1/2 + 2 kappa / beta          (only for long smooth factors)
def worstA(kappa, betas, sigmas, gmax=60, use_pointwise=True):
    b = betas[:, None]; s = sigmas[None, :]
    best = np.ones((len(betas), len(sigmas)))
    for g in range(1, gmax + 1):
        e = np.maximum(1 + g * b * (1 - 2 * s), 2 * g * b * (1 - s))
        best = np.minimum(best, e)
    A = best / (1 - s)
    if use_pointwise:
        A = np.where(s > 0.5 + 2 * kappa / b + 1e-12, 0.0, A)
    i, j = np.unravel_index(np.argmax(A), A.shape)
    return A[i, j], betas[i], sigmas[j]

betas = np.concatenate([np.linspace(0.02, 3.0, 1491), np.linspace(3.02, 8.0, 250)])
sigmas = np.linspace(0.5005, 0.9995, 667)
for kappa in [1/8, 1/6, 0.2, 0.24, 0.249, 0.25]:
    A, bw, sw = worstA(kappa, betas, sigmas)
    target = 2 + 8 * kappa
    # grid sup must not exceed the analytic value (up to grid resolution) and must come close to it
    ok = (A <= target + 2e-3) and (A >= target - 0.03)
    check(f"A(kappa={kappa:.4f}) = sup_grid {A:.4f} vs 2+8kappa = {target:.4f}", ok,
          f"worst beta={bw:.3f} (theory 1+4kappa={1+4*kappa:.3f}), sigma={sw:.4f} (theory {0.5+2*kappa/(1+4*kappa):.4f})")
# the long-polynomial sup restricted to beta in [1/2+2kappa, 1) equals 2+8kappa at the left end (formula f)
k, bb = sp.symbols('kappa beta', positive=True)
f = 2 * bb + 2 * bb * (1 - bb) / (bb - 4 * k)
check("identity f(beta) = 2-8kappa+8kappa(1-4kappa)/(beta-4kappa)",
      sp.simplify(f - (2 - 8 * k + 8 * k * (1 - 4 * k) / (bb - 4 * k))) == 0)
check("f(1/2+2kappa) = 2+8kappa", sp.simplify(f.subs(bb, sp.Rational(1, 2) + 2 * k) - (2 + 8 * k)) == 0)
check("2beta/(beta-4kappa) = 2beta at beta = 1+4kappa (value 2+8kappa)",
      sp.simplify((2 * bb / (bb - 4 * k) - 2 * bb).subs(bb, 1 + 4 * k)) == 0
      and sp.simplify((2 * bb).subs(bb, 1 + 4 * k) - (2 + 8 * k)) == 0)
# short polynomials (no pointwise bound): N^g just above M gives exponent 2 g beta(1-s) <= 2(1+beta)(1-s)
A_short, _, _ = worstA(0.25, np.linspace(0.005, 0.05, 46), sigmas, gmax=400, use_pointwise=False)
check("short factors (beta <= 0.05): sup A <= 2(1+0.05) = 2.1", A_short <= 2.1 + 1e-9, f"{A_short:.4f}")
# without the pointwise bound the method gives no exponent below 4 (pointwise bound is essential)
A_np, bnp, snp = worstA(0.25, betas, sigmas, use_pointwise=False)
check("without the pointwise bound the sup is not below 4 (it grows like 2*beta_max: pointwise bound essential)", A_np >= 3.99, f"{A_np:.3f} at beta={bnp:.3f}")
# exact values
for kap, th in [(sp.Rational(1, 6), sp.Rational(3, 10)), (sp.Rational(1, 4), sp.Rational(1, 4))]:
    check(f"theta(kappa={kap}) = 1/(2+8kappa) = {th}", sp.Rational(1, 1) / (2 + 8 * kap) == th)
check("theta(kappa) > 1/4  <=>  kappa < 1/4 (kappa = 1/4 - delta)",
      sp.simplify(1 / (2 + 8 * (sp.Rational(1, 4) - sp.Symbol('d', positive=True))) - 1 / (4 - 8 * sp.Symbol('d', positive=True))) == 0)
# comparison with the referee's Stucky-set scan at L = 1 (values in referee/kappa_scan_out.txt)
for kappa, refA in [(1/6, 3.3327), (0.24, 3.9182), (0.25, 3.9983)]:
    A, _, _ = worstA(kappa, betas, sigmas)
    check(f"agrees with referee scan at L=1, kappa={kappa:.4f}: {A:.4f} vs {refA}", abs(A - refA) < 0.02)

# ---------------------------------------------------------------------------------------------
print("== B. mean-value lemma C' (Lemma 5)")
def ideals_upto(N):
    out = []
    r = int(math.isqrt(N))
    for a in range(1, r + 1):
        for b in range(0, r + 1):
            n = a * a + b * b
            if n <= N:
                out.append((a, b))
    return out

def gfactor_mu_d(a, b):
    """Moebius and divisor function d_2 of the ideal (a+bi) in Z[i], via the norm factorisation."""
    n = a * a + b * b
    g = math.gcd(a, b)
    fn = sp.factorint(n)
    fg = sp.factorint(g)
    mu, d = 1, 1
    for p, e in fn.items():
        v = fg.get(p, 0)           # v_p(content)
        if p == 2:
            ex = e                 # (1+i)^e
            mu *= 0 if ex >= 2 else -1
            d *= ex + 1
        elif p % 4 == 3:
            ex = e // 2            # (p)^(e/2)
            mu *= 0 if ex >= 2 else -1
            d *= ex + 1
        else:
            # exponents of pi and pibar: v + w and v, where w = e - 2v sits on one side
            e1, e2 = e - v, v
            for ex in (e1, e2):
                if ex >= 2:
                    mu = 0
                elif ex == 1:
                    mu = -mu
            d *= (e1 + 1) * (e2 + 1)
    return mu, d

def mv_test(N, coeffs_fn, M, T, label):
    I = ideals_upto(N)
    data = []
    for (a, b) in I:
        c = coeffs_fn(a, b)
        if c != 0:
            data.append((a, b, c))
    if not data:
        return None
    A_ = np.array([t[0] for t in data], float); B_ = np.array([t[1] for t in data], float)
    c = np.array([t[2] for t in data], complex)
    nrm = A_ ** 2 + B_ ** 2
    x = (2 / math.pi) * np.arctan2(B_, A_)          # lambda^m(a) = e(m x_a)
    u = np.log(nrm)                                  # Na^{-it} = e^{-i t u}
    g = np.array([math.gcd(int(p), int(q)) for p, q in zip(A_, B_)], float)
    Nsup = nrm.max()
    D = np.sqrt(Nsup / (nrm / g ** 2))
    lhs = 0.0
    n = len(c)
    for i0 in range(0, n, 400):
        xi = x[i0:i0 + 400, None] - x[None, :]
        ui = u[i0:i0 + 400, None] - u[None, :]
        s = np.sin(math.pi * xi)
        with np.errstate(divide='ignore', invalid='ignore'):
            K = np.where(np.abs(s) < 1e-13, 2 * M + 1, np.sin((2 * M + 1) * math.pi * xi) / s)
            J = np.where(np.abs(ui) < 1e-13, 2 * T, 2 * np.sin(T * ui) / ui)
        lhs += np.real(np.sum((c[i0:i0 + 400, None] * np.conj(c[None, :])) * K * J))
    sumc2 = float(np.sum(np.abs(c) ** 2))
    bound = E ** 2 * SQPI * (SQPI * M + math.pi * Nsup / 2) * float(np.sum(np.abs(c) ** 2 * (T + SQPI * D)))
    claimC = (M * T + Nsup * min(M, T)) * sumc2           # 'Ricci C' as stated, arbitrary coefficients
    natural = (M + Nsup) * T * sumc2
    duke = float(np.sum((nrm + M * T) * np.abs(c) ** 2))   # Duke 1989 p.205 conjectured form, q = 1
    mv_test.last_duke = duke
    return lhs, bound, claimC, natural, sumc2, n

random.seed(1)
cache_md = {}
def md(a, b):
    if (a, b) not in cache_md:
        cache_md[(a, b)] = gfactor_mu_d(a, b)
    return cache_md[(a, b)]

N0 = 2500
seqs = {
    "rational integers (counterexample)": lambda a, b: 1.0 if b == 0 else 0.0,
    "Moebius on (N/2,N]": lambda a, b: float(md(a, b)[0]) if a * a + b * b > N0 // 2 else 0.0,
    "divisor d_2 on (N/2,N]": lambda a, b: float(md(a, b)[1]) if a * a + b * b > N0 // 2 else 0.0,
    "1 on all ideals (N/2,N]": lambda a, b: 1.0 if a * a + b * b > N0 // 2 else 0.0,
    "random phases, primitive only": lambda a, b: cmath.exp(2j * math.pi * random.random()) if math.gcd(a, b) == 1 else 0.0,
    "random phases, all ideals": lambda a, b: cmath.exp(2j * math.pi * random.random()),
}
nat_ratio = {}
for label, fn in seqs.items():
    for (M, T) in [(10, 1), (300, 1), (5000, 1), (5000, 4), (100000, 1)]:
        lhs, bound, claimC, natural, sumc2, n = mv_test(N0, fn, M, T, label)
        check(f"C' holds: {label}, N={N0}, M={M}, T={T}: LHS/bound = {lhs/bound:.3e}", lhs <= bound and lhs >= -1e-6 * bound)
        if label.startswith("rational") and M >= 5000:
            check(f"  'Ricci C' as stated fails for the counterexample (M={M},T={T}): LHS/[(MT+N min(M,T))sum|c|^2] = {lhs/claimC:.1f} > 1", lhs > claimC)
        if M == 100000:
            nat_ratio[label] = lhs / natural
            if label.startswith("rational"):
                check(f"  Duke's conjectured bound sum (N(a)+qMT)|c(a)|^2 fails for the counterexample: LHS/RHS = {lhs/mv_test.last_duke:.1f}", lhs > 10 * mv_test.last_duke)
rm = nat_ratio["Moebius on (N/2,N]"]; rc = nat_ratio["rational integers (counterexample)"]
check(f"at M=1e5,T=1: LHS/[(M+N)T sum|c|^2] is {rm:.2f} for Moebius (bounded) but {rc:.1f} for the counterexample (the content weight is needed)", rm <= 8 and rc >= 5 * rm)
# primitive ideals have 1/N-spaced angles mod pi/2: minimal spacing of x_a over primitive ideals
I = [(a, b) for (a, b) in ideals_upto(4000) if math.gcd(a, b) == 1]
xs = np.sort((2 / math.pi) * np.arctan2([q for p, q in I], [p for p, q in I]))
gaps = np.diff(np.concatenate([xs, [xs[0] + 1]]))
check(f"primitive ideals, N<=4000: min gap of x_a = {gaps.min():.3e} >= 2/(pi N) = {2/(math.pi*4000):.3e}", gaps.min() >= 2 / (math.pi * 4000))
# the collapse identity sum_a c(a) lambda^m(a) Na^{-it} = sum_{a'} lambda^m(a') Na'^{-it} C_{a'}(t)
I = ideals_upto(600)
cc = {ab: complex(random.gauss(0, 1), random.gauss(0, 1)) for ab in I}
for (m, t) in [(3, 0.7), (-17, 2.3), (41, -5.1)]:
    lhs = sum(c * cmath.exp(4j * m * math.atan2(b, a)) * (a * a + b * b) ** (-1j * t) for (a, b), c in cc.items())
    rhs = 0
    for (a, b) in I:
        if math.gcd(a, b) != 1:
            continue
        Cap = 0
        d = 1
        while d * d * (a * a + b * b) <= 600:
            Cap += cc[(d * a, d * b)] * d ** (-2j * t)
            d += 1
        rhs += cmath.exp(4j * m * math.atan2(b, a)) * (a * a + b * b) ** (-1j * t) * Cap
    check(f"collapse identity at m={m}, t={t}", abs(lhs - rhs) < 1e-8 * (1 + abs(lhs)))

# ---------------------------------------------------------------------------------------------
print("== C. divisor moments with the content weight (Lemma 6)")
def spf_sieve(n):
    s = np.zeros(n + 1, dtype=np.int32)
    s[1] = 1
    for p in range(2, n + 1):
        if s[p] == 0:
            s[p::p][s[p::p] == 0] = p
    return s
def fact(n, spf):
    out = {}
    while n > 1:
        p = int(spf[n]); e = 0
        while n % p == 0:
            n //= p; e += 1
        out[p] = e
    return out
Y = 60000
spf = spf_sieve(Y)
def d2_ideal(a, b):
    n = a * a + b * b; g = math.gcd(a, b)
    fn = fact(n, spf); fg = fact(g, spf) if g > 1 else {}
    d = 1
    for p, e in fn.items():
        v = fg.get(p, 0)
        if p == 2:
            d *= e + 1
        elif p % 4 == 3:
            d *= e // 2 + 1
        else:
            d *= (e - v + 1) * (v + 1)
    return d
vals = {}
r = int(math.isqrt(Y))
S_w = {}; S_0 = {}
for Yc in [2000, 8000, 30000, 60000]:
    S_w[Yc] = 0.0; S_0[Yc] = 0.0
for a in range(1, r + 1):
    for b in range(0, r + 1):
        n = a * a + b * b
        if n > Y:
            break
        d = d2_ideal(a, b); g = math.gcd(a, b)
        for Yc in S_w:
            if n <= Yc:
                S_w[Yc] += d * d * g / n
                S_0[Yc] += d * d / n
for Yc in sorted(S_w):
    L = math.log(Yc)
    print(f"   Y={Yc}: sum d^2 cont/N = {S_w[Yc]:.2f}, sum d^2/N = {S_0[Yc]:.2f}, ratio {S_w[Yc]/S_0[Yc]:.3f}, sum d^2 cont/N / (log Y)^5 = {S_w[Yc]/L**5:.4f}")
rats = [S_w[Yc] / S_0[Yc] for Yc in sorted(S_w)]
check("content weight costs at most a log factor: ratio(Y) <= 2 log Y for all Y tested",
      all(rt <= 2 * math.log(Yc) for rt, Yc in zip(rats, sorted(S_w))), ", ".join(f"{x:.2f}" for x in rats))
# product bound used in the proof: sum_{Na<=Y} d(a)^2 cont(a)/Na <= (sum_{d<=sqrt Y} d((d))^2/d)(sum_{Na<=Y} d(a)^2/Na)
def d2_rational_ideal(dd):
    if dd == 1:
        return 1
    out = 1
    for p, e in fact(dd, spf).items():
        out *= (2 * e + 1) if p == 2 else ((e + 1) if p % 4 == 3 else (e + 1) ** 2)
    return out
for Yc in sorted(S_w):
    S1 = sum(d2_rational_ideal(dd) ** 2 / dd for dd in range(1, int(math.isqrt(Yc)) + 1))
    check(f"product bound at Y={Yc}: {S_w[Yc]:.1f} <= {S1:.2f} * {S_0[Yc]:.1f}", S_w[Yc] <= S1 * S_0[Yc])
# submultiplicativity d_g(ab) <= d_g(a) d_g(b) for prime powers (binomial inequality), g <= 6, exponents <= 12
okb = all(math.comb(a + b + g - 1, g - 1) <= math.comb(a + g - 1, g - 1) * math.comb(b + g - 1, g - 1)
          for g in range(1, 7) for a in range(13) for b in range(13))
check("d_g submultiplicative on prime powers (g<=6, exponents<=12)", okb)

# ---------------------------------------------------------------------------------------------
print("== D. lattice points in polar boxes (Lemma 2)")
random.seed(7)
okD = True; worst = 0
for _ in range(150):
    r1 = random.uniform(2, 250); y = random.uniform(1, 3 * r1 * r1)
    r2 = math.sqrt(r1 * r1 + y)
    om = random.uniform(1e-4, 0.6); ph = random.uniform(-math.pi, math.pi)
    R = int(r2) + 1
    a = np.arange(-R, R + 1)[:, None]; b = np.arange(-R, R + 1)[None, :]
    rr = np.sqrt(a * a + b * b)
    ang = np.mod(np.arctan2(b, a) - ph, 2 * math.pi)
    cnt = int(np.sum((rr >= r1) & (rr <= r2) & (ang <= om)))
    bnd = (om / 2 + 1.12 / r1) * (y + 2 * (r1 + r2))
    worst = max(worst, cnt / bnd)
    okD &= cnt <= bnd
check(f"polar box count <= (omega/2 + 1.12/r1)(r2^2-r1^2 + 2(r1+r2)) on 150 random boxes (max ratio {worst:.3f})", okD)
check("arcsin(1/(sqrt2 r)) <= 1.111/r for r >= 2", all(math.asin(1 / (math.sqrt(2) * r)) <= 1.111 / r for r in np.linspace(2, 1e4, 20001)))

# ---------------------------------------------------------------------------------------------
print("== E. angular cut-off (Lemma 1)")
alpha, beta, Dl, rr_ = 0.20, 0.23, 0.006, 3
def ghat(kk):
    if kk == 0:
        return beta - alpha
    w = beta - alpha
    return (math.sin(math.pi * kk * w) / (math.pi * kk)) * (math.sin(math.pi * kk * Dl / rr_) / (math.pi * kk * Dl / rr_)) ** rr_
okE = True
for kk in list(range(1, 2000)) + [10 ** 4, 10 ** 5, 10 ** 6]:
    bnd = min(beta - alpha, 1 / (math.pi * kk)) * min(1.0, (rr_ / (math.pi * kk * Dl)) ** rr_)
    okE &= abs(ghat(kk)) <= bnd * (1 + 1e-12)
check("|g_k| <= min(beta-alpha, 1/(pi|k|)) min(1, (r/(pi|k|Delta))^r) for k <= 2000 and 1e4..1e6", okE)
# direct convolution on a grid: plateau, support, range
n = 2 ** 18
xg = np.arange(n) / n
ind = ((xg >= alpha) & (xg <= beta)).astype(float)
box_w = Dl / rr_
box = ((xg < box_w / 2) | (xg > 1 - box_w / 2)).astype(float)
box /= box.sum()
psi = ind.copy()
for _ in range(rr_):
    psi = np.real(np.fft.ifft(np.fft.fft(psi) * np.fft.fft(box)))
tol = 3.0 / (n * box_w)
plateau = (xg >= alpha + Dl / 2 + 2 / n) & (xg <= beta - Dl / 2 - 2 / n)
outside = (xg <= alpha - Dl / 2 - 2 / n) | (xg >= beta + Dl / 2 + 2 / n)
check("psi = 1 on [alpha+Delta/2, beta-Delta/2], 0 outside [alpha-Delta/2, beta+Delta/2], 0<=psi<=1",
      np.all(np.abs(psi[plateau] - 1) < tol) and np.all(np.abs(psi[outside]) < tol) and psi.min() > -tol and psi.max() < 1 + tol)
F = np.fft.fft(psi) / n
okF = all(abs(abs(F[kk]) - abs(ghat(kk))) < 2e-4 for kk in range(0, 60))
check("Fourier coefficients of the convolved function match the closed form (|.|, k<60)", okF)

# ---------------------------------------------------------------------------------------------
print("== F. Gaussian norm weight (Lemma 3)")
mp.mp.dps = 30
x0, H0 = mp.mpf(100), mp.mpf(3)
Vf = lambda u: mp.erfc(H0 * mp.log(u / x0) / mp.sqrt(2)) / 2
for s in [mp.mpc(0.5, 2), mp.mpc(1.2, -0.7), mp.mpc(0.5, 9)]:
    num = mp.quad(lambda v: Vf(x0 * mp.e ** v) * (x0 * mp.e ** v) ** s, [-mp.inf] + [mp.mpf(t) for t in np.linspace(-60, 4, 129)] + [mp.inf])
    clo = x0 ** s * mp.e ** (s * s / (2 * H0 ** 2)) / s
    check(f"Mellin transform of V at s={complex(s)}: rel.err {float(abs(num-clo)/abs(clo)):.1e}", abs(num - clo) / abs(clo) < 1e-10)
okT = all(float(mp.erfc(yv / mp.sqrt(2)) / 2) <= math.exp(-yv * yv / 2) for yv in np.linspace(0, 12, 241))
check("P(Z >= y) <= exp(-y^2/2) for y in [0,12]", okT)

# ---------------------------------------------------------------------------------------------
print("== G. Heath-Brown identity on ideals of Z[i]")
def canon(a, b):
    for _ in range(4):
        if a > 0 and b >= 0:
            return (a, b)
        a, b = -b, a
    raise ValueError
Ymax = 220
IDs = sorted(ideals_upto(Ymax), key=lambda t: (t[0] ** 2 + t[1] ** 2, t))
Nm = {ab: ab[0] ** 2 + ab[1] ** 2 for ab in IDs}
def mul(p, q):
    return canon(p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0])
pairs = [(p, q) for p in IDs for q in IDs if Nm[p] * Nm[q] <= Ymax]
def conv(f, g):
    h = {ab: 0.0 for ab in IDs}
    for p, q in pairs:
        fp = f.get(p, 0.0)
        if fp == 0.0:
            continue
        gq = g.get(q, 0.0)
        if gq != 0.0:
            h[mul(p, q)] += fp * gq
    return h
MU, LAM = {}, {}
for ab in IDs:
    mu, _ = gfactor_mu_d(*ab) if ab != (1, 0) else (1, 1)
    MU[ab] = float(mu)
    n = Nm[ab]
    fa = sp.factorint(n) if n > 1 else {}
    # Lambda of the ideal: log N(prime ideal) if the ideal is a prime power
    lam = 0.0
    if n > 1:
        if len(fa) == 1:
            p, e = list(fa.items())[0]
            g = math.gcd(*ab)
            if p == 2:
                lam = math.log(2)
            elif p % 4 == 3:
                lam = math.log(p * p) if (e % 2 == 0 and (ab[1] == 0)) else 0.0
            else:
                lam = math.log(p) if g == 1 else 0.0
    LAM[ab] = lam
ONE = {ab: 1.0 for ab in IDs}
LOG = {ab: math.log(Nm[ab]) for ab in IDs}
for (X, J) in [(14, 2), (6, 3)]:
    rng = random.Random(J)
    for smooth in [False, True]:
        wX = {}
        for ab in IDs:
            t = Nm[ab] / X
            wX[ab] = 1.0 if t <= 1 else ((rng.random() if t < 2 else 0.0) if smooth else 0.0)
        muX = {ab: MU[ab] * wX[ab] for ab in IDs}
        total = {ab: 0.0 for ab in IDs}
        for j in range(1, J + 1):
            h = dict(LOG)
            for _ in range(j - 1):
                h = conv(h, ONE)
            for _ in range(j):
                h = conv(h, muX)
            for ab in IDs:
                total[ab] += (-1) ** (j - 1) * math.comb(J, j) * h[ab]
        err = max(abs(total[ab] - LAM[ab]) for ab in IDs if Nm[ab] <= X ** J)
        check(f"HB identity, X={X}, J={J}, {'smoothed' if smooth else 'sharp'} Moebius cut: max error on N a <= X^J = {err:.1e}", err < 1e-9)
# Lambda sanity: sum_{Na <= 220} Lambda(a) vs psi_K
check("sum of Lambda over ideals of norm <= 220 is ~ 220 (prime ideal theorem scale)", 150 < sum(LAM.values()) < 300, f"{sum(LAM.values()):.1f}")

# ---------------------------------------------------------------------------------------------
print("== H. characters mod q, units, infinity types")
def residues_units(qa, qb):
    """Units of Z[i]/(q), q = qa + qb i: classes keyed by z*conj(q) mod N(q) (a complete invariant)."""
    Nq = qa * qa + qb * qb
    def key(a, b):
        re_, im_ = a * qa + b * qb, b * qa - a * qb          # (a+bi)(qa - qb i)
        return (re_ % Nq, im_ % Nq)
    reps = {}
    for a in range(Nq):
        for b in range(Nq):
            reps.setdefault(key(a, b), (a, b))
    assert len(reps) == Nq
    units = [z for z in reps.values()
             if any(key(z[0] * w[0] - z[1] * w[1] - 1, z[0] * w[1] + z[1] * w[0]) == (0, 0) for w in reps.values())]
    return units
for (qa, qb, Phi_expected) in [(2, 1, 4), (5, 0, 16), (3, 0, 8)]:
    U = residues_units(qa, qb)
    # residues coprime to q have norm prime to N(q) only when q is a product of split primes or inert; check count
    check(f"Phi(q) for q={qa}+{qb}i: {len(U)} (expected {Phi_expected})", len(U) == Phi_expected)
# the sum over units sum_eps chi(eps) eps^k = 4 [chi(i) i^k = 1]
for chi_i in [1, 1j, -1, -1j]:
    for kk in range(-8, 9):
        sm = sum((chi_i ** e) * (1j ** e) ** kk for e in range(4))
        expect = 4 if abs(chi_i * 1j ** kk - 1) < 1e-12 else 0
        if abs(sm - expect) > 1e-9:
            check("unit sum", False)
check("sum over units = 4*[chi(i) i^k = 1] for all chi(i) in mu_4 and |k| <= 8", True)
# MV conductor arithmetic: C(xi |.|^{it}) <= N(q) (2 + |k|/2 + |t|)^2 and C^{1/4-delta} <= c (1+|k|+|t|)^{1/2-2delta}
dlt = 0.01
okC = all((5 * (2 + kk / 2 + tt) ** 2) ** (0.25 - dlt) <= 5 ** 0.25 * 2 * (1 + kk + tt) ** (0.5 - 2 * dlt)
          for kk in range(0, 2000, 7) for tt in np.linspace(0, 50, 11))
check("conductor bound (N(q)(2+|k|/2+|t|)^2)^(1/4-d) <= 2 N(q)^(1/4) (1+|k|+|t|)^(1/2-2d)", okC)

# ---------------------------------------------------------------------------------------------
print("== I. Gaussian primes in residue classes in sectors of width c x^(-1/4), x = 1e7 (sanity check of the main term)")
X = 10 ** 7
isp = np.ones(X + 1, dtype=bool); isp[:2] = False
for p in range(2, int(X ** 0.5) + 1):
    if isp[p]:
        isp[p * p::p] = False
base = np.zeros(X + 1, dtype=np.int32)
pr = np.nonzero(isp)[0]
base[pr] = pr
for p in pr[pr <= int(X ** 0.5) + 1]:
    q = int(p) * int(p)
    while q <= X:
        base[q] = p
        q *= int(p)
del pr
def lam_points(a, b):
    n = a * a + b * b
    p = base[n]
    lam = np.zeros(len(n))
    pos = p > 0
    prime = pos & (p == n)
    lam[prime] = np.log(n[prime].astype(float))
    pw = pos & (p != n)
    g = np.gcd(a, b)
    two = pw & (p == 2)
    lam[two] = math.log(2)
    one = pw & (p % 4 == 1) & (g == 1)
    lam[one] = np.log(p[one].astype(float))
    three = pw & (p % 4 == 3) & ((a == 0) | (b == 0))
    lam[three] = 2 * np.log(p[three].astype(float))
    return lam
def sector_points(phi, om, X):
    R = int(math.isqrt(X))
    As, Bs = [], []
    lo = phi
    for a in range(-R, R + 1):
        bmax = int(math.isqrt(X - a * a))
        b = np.arange(-bmax, bmax + 1)
        ang = np.mod(np.arctan2(b, a) - lo, 2 * math.pi)
        sel = (ang > 0) & (ang <= om)
        if sel.any():
            As.append(np.full(sel.sum(), a)); Bs.append(b[sel])
    return np.concatenate(As), np.concatenate(Bs)
zs = []
pooled = []
for (phi, c, orient) in [(0.31, 1.0, +1), (0.77, 1.0, -1), (1.13, 2.0, +1), (2.40, 1.0, +1), (-1.95, 1.5, -1), (3.0, 4.0, +1)]:
    om = c * X ** (-0.25)
    lo = phi if orient > 0 else phi - om
    a, b = sector_points(lo, om, X)
    lam = lam_points(a, b)
    tot = lam.sum()
    main_all = (2 / math.pi) * om * X
    pooled.append((tot - main_all) / math.sqrt((lam ** 2).sum()))
    # q = 2+i: class = (a - 2b) mod 5, Phi = 4
    cls = np.mod(a - 2 * b, 5)
    for r_ in range(1, 5):
        s = lam[cls == r_].sum(); v = (lam[cls == r_] ** 2).sum()
        zs.append((s - main_all / 4) / math.sqrt(v))
    # q = 5: class = (a mod 5, b mod 5) with N coprime to 5, Phi = 16
    nn = (a * a + b * b) % 5
    key = np.mod(a, 5) * 5 + np.mod(b, 5)
    for kq in range(25):
        if (kq // 5) ** 2 + (kq % 5) ** 2 == 0:
            continue
        sel = (key == kq) & (nn != 0)
        if not sel.any():
            continue
        s = lam[sel].sum(); v = (lam[sel] ** 2).sum()
        zs.append((s - main_all / 16) / math.sqrt(v))
    print(f"   sector phi={phi:+.2f} width {c}x^(-1/4) orient {orient:+d}: points {len(a)}, sum Lambda {tot:.0f}, main (2/pi) omega x = {main_all:.0f}, ratio {tot/main_all:.4f}")
zs = np.array(zs)
print(f"   {len(zs)} (sector, class) z-scores: mean {zs.mean():+.3f}, rms {np.sqrt((zs**2).mean()):.3f}, max |z| {np.abs(zs).max():.2f}")
check("class sums mod 2+i and mod 5 agree with (2/pi) omega x / Phi(q): rms z in [0.5, 1.6], max |z| < 4.5",
      0.5 <= np.sqrt((zs ** 2).mean()) <= 1.6 and np.abs(zs).max() < 4.5)
check("pooled sector sums agree with (2/pi) omega x: max |z| < 4", max(abs(p) for p in pooled) < 4, ", ".join(f"{p:+.2f}" for p in pooled))
del base, isp

# ---------------------------------------------------------------------------------------------
print("== J. parameter bookkeeping")
kappa, theta, eps, Dlt, J = 0.24, 0.252, 1e-4, 1e-3, 20
A = 2 + 8 * kappa
c1 = (theta + 2 * eps) * A < 1 - 3 * Dlt
c2 = 2 * (theta + 2 * eps + 1 / J) < 1 - 3 * Dlt
c3 = theta + 2 * eps + 2 * Dlt < 0.5
check(f"kappa=0.24: theta=0.252 admissible: (theta+2eps)A={(theta+2*eps)*A:.4f}<1-3Delta; short: {2*(theta+2*eps+1/J):.4f}<1-3Delta; I_0 case", c1 and c2 and c3)
check("theta = 1/4 is admissible for every kappa < 1/4 (strict inequality with room)", all((0.25) * (2 + 8 * kp) < 1 for kp in np.linspace(0.125, 0.2499, 100)))

# ---------------------------------------------------------------------------------------------
print("== K. large sieve (Gaussian majorant form) and Gallagher's inequality")
rng = np.random.default_rng(3)
okK = True
for trial in range(20):
    R = 200
    xs = np.sort(rng.random(R))
    dlt = np.min(np.diff(np.concatenate([xs, [xs[0] + 1]])))
    Mm = int(rng.integers(5, 400))
    bvec = rng.normal(size=R) + 1j * rng.normal(size=R)
    ms = np.arange(-Mm, Mm + 1)
    S = np.exp(2j * math.pi * np.outer(ms, xs)) @ bvec
    lhs = np.sum(np.abs(S) ** 2)
    okK &= lhs <= E * (SQPI * Mm + 1 / dlt) * np.sum(np.abs(bvec) ** 2)
check("sum_{|m|<=M} |sum_r b_r e(m x_r)|^2 <= e(sqrt(pi) M + 1/delta) sum|b_r|^2 (20 random cases)", okK)
okG = True
for trial in range(20):
    n_ = rng.integers(2, 60, size=30); cf = rng.normal(size=30) + 1j * rng.normal(size=30)
    t0 = rng.uniform(-50, 50)
    ts = np.linspace(t0 - 0.5, t0 + 0.5, 4001)
    Pt = np.exp(-1j * np.outer(ts, np.log(n_))) @ cf
    dPt = np.exp(-1j * np.outer(ts, np.log(n_))) @ (cf * (-1j) * np.log(n_))
    rhs = np.trapezoid(np.abs(Pt) ** 2 + 2 * np.abs(Pt) * np.abs(dPt), ts)
    P0 = np.exp(-1j * t0 * np.log(n_)) @ cf
    okG &= abs(P0) ** 2 <= rhs * (1 + 1e-6)
check("Gallagher: |P(t0)|^2 <= int_{|t-t0|<=1/2} (|P|^2 + 2|P||P'|) (20 random Dirichlet polynomials)", okG)

print(f"{NPASS + NFAIL} checks, {NFAIL} failed")
