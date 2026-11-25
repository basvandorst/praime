"""check_arcs.py: numerical checks for arcs_theorem.md (major arcs of the radial resonance of Gaussian primes).

Sections (each prints PASS/FAIL lines; the script exits with an AssertionError on the first FAIL):
  A. algebra (Lemma 2): for w primitive, |w| = h, zeta -> Re(zeta conj w) mod h maps (Z[i]/h)^* onto (Z/h)^* with
     equal fibres; mean of e(-r/h) = mu(h)/phi(h); count of primitive w with |w| = h is 4 * 2^omega(h).
  B. identity (Corollary 1): sum over rays of A_h (w/h)^m (normalised per axis) = prod_p (1 - 2cos(2m theta_p)/(p-1)),
     term by term in h (ray enumeration vs multiplicative formula), partial sums to h <= 10^4 vs the product to 10^7.
  C. Fresnel factor Phi_M = int F(y/M) e^{-i pi y^2} dy vs e^{-i pi/4}; closed form W_V Phi_M X^{3/4} vs quadrature.
  D. zone sums over Gaussian primes with N z <= 10^7 (m = 4 and m = 0, M = 3) vs the theorem's main term
     A_h u^m int G dA, plus a heuristic overlap model (neighbouring rays), noise from the random-phase model and from
     four control frequencies.
  E. elementary statement: Re(z conj w) mod h for Gaussian primes in the flat part of each zone (never 0 mod p | h;
     uniform on (Z/h)^*; Lambda-weighted mean of e(-r/h) vs mu(h)/phi(h)); lattice-point control for 3 + 4i.
Run: (ulimit -v 1200000; timeout 600 nice python -I -B check_arcs.py)
"""
import math, cmath, sys, time
from math import gcd
import numpy as np

T0 = time.time()
TWO_PI = 2 * math.pi
X = 10 ** 7          # norms N z <= X
M = 3.0              # zone parameter
m = 4                # Hecke character psi_m(z) = (z/|z|)^m
V_LO, V_HI, V_W = 0.05, 1.0, 0.25    # V supported on [0.05, 1] (units of X), flat on [0.30, 0.75]
CONTROLS = (0.6713, 0.8219, 1.1597, 1.3711)
FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ------------------------------------------------------------------ smooth cutoffs
def _f(x):
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    pos = x > 0
    out[pos] = np.exp(-1.0 / x[pos])
    return out


def step(x):
    """C-infinity step: 0 for x <= 0, 1 for x >= 1."""
    a = _f(x); b = _f(1.0 - np.asarray(x, float))
    return a / (a + b)


def F(y):
    """even, 1 on [-1, 1], 0 outside [-2, 2]."""
    return step(2.0 - np.abs(y))


def V(v):
    return step((np.asarray(v, float) - V_LO) / V_W) * step((V_HI - np.asarray(v, float)) / V_W)


def ee(x):
    return np.exp(1j * TWO_PI * x)


# ------------------------------------------------------------------ small number theory
def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def mobius(n):
    f = factor(n)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)


def phi(n):
    r = n
    for p in factor(n):
        r = r // p * (p - 1)
    return r


def prim_w(h):
    """all primitive w = (x, y) in Z^2 with x^2 + y^2 = h^2 (gcd(x, y) = 1)."""
    out = []
    for x in range(-h, h + 1):
        y2 = h * h - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            for yy in ({y, -y} if y else {0}):
                if gcd(x, yy) == 1:
                    out.append((x, yy))
    return out


# ================================================================== A. algebra
print("== A. algebra (Lemma 2)")
for h in (1, 5, 13, 17, 25, 29, 37, 65, 85, 125, 145):
    ws = prim_w(h)
    om = len(factor(h)) if h > 1 else 0
    check(len(ws) == 4 * 2 ** om, f"h={h}: {len(ws)} primitive w with |w|=h (expected 4*2^omega = {4 * 2 ** om})")
    units = [(a, b) for a in range(h) for b in range(h)
             if all((a * a + b * b) % p for p in factor(h))] if h > 1 else [(0, 0)]
    Phi = len(units)
    for (wx, wy) in ws:
        cnt = {}
        for (a, b) in units:
            r = (a * wx + b * wy) % h
            cnt[r] = cnt.get(r, 0) + 1
        img_ok = all(gcd(r, h) == 1 for r in cnt) and len(cnt) == phi(h)
        fib_ok = len(set(cnt.values())) == 1
        mean = sum(c * cmath.exp(-2j * math.pi * r / h) for r, c in cnt.items()) / Phi
        ok = img_ok and fib_ok and abs(mean - mobius(h) / phi(h)) < 1e-12
        if not ok:
            check(False, f"h={h}, w={wx}+{wy}i: image/fibres/mean wrong (mean {mean})")
            break
    else:
        check(True, f"h={h}: all {len(ws)} w: image (Z/h)^*, equal fibres {Phi // phi(h)}, mean e(-r/h) = mu/phi = {mobius(h) / phi(h):+.6f}")

# ================================================================== B. identity
print("== B. identity (Corollary 1), m = %d" % m)
HMAX = 10 ** 4
R = math.isqrt(HMAX) + 1
rays = {}            # h -> list of complex (w/h)^m for all primitive w with |w| = h
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        n = a * a + b * b
        if n == 0 or n > HMAX or n % 2 == 0 or gcd(a, b) != 1:
            continue
        w = complex(a, b) ** 2
        for e in (1, 1j, -1, -1j):
            rays.setdefault(n, set()).add((round((e * w).real), round((e * w).imag)))
rays[1] = {(1, 0), (0, 1), (-1, 0), (0, -1)}
term_rays = {}
for h, ws in rays.items():
    s = sum(complex(x, y) ** m for (x, y) in ws) / h ** m
    term_rays[h] = (mobius(h) / phi(h)) * s / 4.0          # A_h/(4/pi) summed over the rays, per axis (4 axes -> 1)
# multiplicative formula: g(h) = prod_{p | h} (-2 cos(2 m theta_p)/(p - 1)) for squarefree h, p = 1 mod 4
theta = {}
for a in range(1, R + 1):
    for b in range(1, a):
        p = a * a + b * b
        if p <= HMAX and len(factor(p)) == 1 and factor(p).get(p) == 1:
            theta[p] = math.atan2(b, a)
maxdiff = 0.0
for h in range(1, HMAX + 1):
    f = factor(h) if h > 1 else {}
    if h > 1 and (any(e > 1 for e in f.values()) or any(p % 4 != 1 for p in f)):
        g = 0.0
    else:
        g = 1.0
        for p in f:
            g *= -2 * math.cos(2 * m * theta[p]) / (p - 1)
    tr = term_rays.get(h, 0.0)
    maxdiff = max(maxdiff, abs(tr - g))
check(maxdiff < 1e-12, f"ray enumeration = multiplicative formula for every h <= {HMAX} (max diff {maxdiff:.1e})")
S = {}
acc = 0.0
for h in range(1, HMAX + 1):
    acc += complex(term_rays.get(h, 0.0)).real
    if h in (10, 100, 1000, 10000):
        S[h] = acc
imag_max = max(abs(complex(v).imag) for v in term_rays.values())
check(imag_max < 1e-12, f"every ray term is real (sum over w and conj w): max imag {imag_max:.1e}")

# ================================================================== C. Fresnel, closed form
print("== C. Fresnel factor and closed form")
yy = np.linspace(-40, 40, 400001)
for MM in (1.0, 2.0, 3.0):
    PhiM = np.trapezoid(F(yy / MM) * np.exp(-1j * math.pi * yy ** 2), yy)
    print(f"   Phi_M (M={MM:.0f}) = {PhiM.real:+.10f} {PhiM.imag:+.10f}i   |Phi_M - e^(-i pi/4)| = {abs(PhiM - cmath.exp(-1j * math.pi / 4)):.2e}")
Phi3 = np.trapezoid(F(yy / M) * np.exp(-1j * math.pi * yy ** 2), yy)
check(abs(Phi3 - cmath.exp(-1j * math.pi / 4)) < 1e-4, "Phi_3 = e^(-i pi/4) to 1e-4 (errors 5e-2, 1e-3, 2e-5 at M = 1, 2, 3: faster than any power)")
vv = np.linspace(0, 1.2, 1200001)
WV = 0.5 * np.trapezoid(np.where(vv > 0, vv, 1) ** -0.25 * V(vv), vv)
print(f"   W_V = (1/2) int v^(-1/4) V(v) dv = {WV:.10f}")


def Jm(phi_u, phi_x, Xv=X, nr=400, dy=0.008):
    """{m: int psi_m(z) e(Re(z conj xi) - |z|) W_u(z) dA} for m = 0, 4; xi = e^{i phi_x}, W_u = V(|z|^2/Xv) F(s_u/(M sqrt t_u)).
    Coordinates z = r e^{i(phi_x + y/sqrt r)}, dA = sqrt(r) dr dy; y over the support of W_u (|y - y0| <= 2.4 M)."""
    r = np.linspace(math.sqrt(V_LO * Xv), math.sqrt(V_HI * Xv), nr)
    tot = {0: 0j, 4: 0j}
    for i0 in range(0, nr - 1, 50):                 # chunks share their end points: no interval is lost
        rc = r[i0:i0 + 51]
        rr = rc[:, None]
        y0 = (phi_u - phi_x) * np.sqrt(rr)
        ny = int(round(2 * 2.4 * M / dy)) + 1
        y = y0 + np.linspace(-2.4 * M, 2.4 * M, ny)[None, :]
        al = y / np.sqrt(rr)
        du = phi_x + al - phi_u
        tu = rr * np.cos(du); su = rr * np.sin(du)
        Wu = V(rr ** 2 / Xv) * np.where(tu > 0, F(su / (M * np.sqrt(np.maximum(tu, 1e-9)))), 0.0)
        base = ee(-2 * rr * np.sin(al / 2) ** 2) * Wu * np.sqrt(rr)
        dyy = y[0, 1] - y[0, 0]
        for mm in (0, 4):
            f = base * np.exp(1j * mm * (phi_x + al)) if mm else base
            tot[mm] += np.trapezoid(np.trapezoid(f, dx=dyy, axis=1), rc)
    return tot


def own_main(mm, Xv=X):
    """int G dA in zone coordinates (t, y = s/sqrt t), G = e^{i m alpha} e(-(|z| - t)) V(N/X) F(y/M):
    an independent quadrature of the theorem's integral."""
    t = np.linspace(math.sqrt(V_LO * Xv) * 0.97, math.sqrt(V_HI * Xv), 1401)
    y = np.linspace(-2 * M, 2 * M, 2401)
    tot = 0j
    for i0 in range(0, len(t) - 1, 100):
        tc = t[i0:i0 + 101]
        tt = tc[:, None]
        s = y[None, :] * np.sqrt(tt)
        nrm = tt ** 2 + s ** 2
        G = np.exp(1j * mm * np.arctan2(s, tt)) * ee(-(np.sqrt(nrm) - tt)) * V(nrm / Xv) * F(y[None, :] / M) * np.sqrt(tt)
        tot += np.trapezoid(np.trapezoid(G, y, axis=1), tc)
    return tot


I_closed = WV * Phi3 * X ** 0.75
I_q4 = own_main(4)
I_q0 = own_main(0)
I_p = Jm(0.0, 0.0)
print(f"   int G dA (m=4): quadrature in (t,y) {I_q4.real:.2f} {I_q4.imag:+.2f}i; in (r,y) {I_p[4].real:.2f} {I_p[4].imag:+.2f}i;"
      f" closed form W_V Phi_M X^3/4 = {I_closed.real:.2f} {I_closed.imag:+.2f}i")
rel4 = abs(I_q4 - I_closed) / abs(I_closed); rel0 = abs(I_q0 - I_closed) / abs(I_closed)
check(rel4 < 2e-3 and rel0 < 2e-3,
      f"int G dA = W_V Phi_M X^(3/4) (1 + O(M^2 X^-1/2)): rel. diff m=4 {rel4:.1e}, m=0 {rel0:.1e}")
check(abs(I_q4 - I_p[4]) / abs(I_q4) < 1e-4 and abs(I_q0 - I_p[0]) / abs(I_q0) < 1e-4,
      f"two quadratures of int G dA agree (rel {abs(I_q4 - I_p[4]) / abs(I_q4):.1e}, {abs(I_q0 - I_p[0]) / abs(I_q0):.1e})")

# ================================================================== D. Gaussian primes up to 10^7
print("== D. Gaussian primes with N z <= %d" % X)
isp = np.ones(X + 1, dtype=bool); isp[:2] = False
for q in range(2, math.isqrt(X) + 1):
    if isp[q]:
        isp[q * q::q] = False
As, Bs = [], []
for a in range(1, math.isqrt(X) + 1):
    Bm = math.isqrt(X - a * a)
    b = np.arange(-Bm, Bm + 1, dtype=np.int64)
    n = a * a + b * b
    sel = isp[n] & (b != 0)
    if isp[a] and a % 4 == 3:
        sel |= (b == 0)
    As.append(np.full(int(sel.sum()), a, np.int64)); Bs.append(b[sel])
A = np.concatenate(As); B = np.concatenate(Bs); del As, Bs
N = A * A + B * B
inert = B == 0
N[inert] = A[inert] ** 2          # (already) norm of the inert prime q is q^2
LAM = np.log(N.astype(float))
RAD = np.sqrt(N.astype(float))
ARG = np.arctan2(B.astype(float), A.astype(float))
print(f"   {len(A)} Gaussian prime elements with Re z > 0 ({time.time() - T0:.0f} s)")
# product prod_{p <= X, p = 1 mod 4} (1 - 2 cos(2 m theta_p)/(p - 1)) from the octant a > b > 0
oc = (B > 0) & (B < A)
pr = N[oc].astype(float); th = np.arctan2(B[oc].astype(float), A[oc].astype(float))
logprod = np.sum(np.log1p(-2 * np.cos(2 * m * th) / (pr - 1)))
PROD = math.exp(logprod)
c_rep = 2.026446079853554     # constants.json, A4, c_exact = B/L(1, Ad)
print(f"   partial singular series: S_10={S[10]:.6f} S_100={S[100]:.6f} S_1000={S[1000]:.6f} S_10000={S[10000]:.6f};"
      f" product p<=1e7: {PROD:.6f}; (pi/4) c_exact = {math.pi / 4 * c_rep:.6f}")
check(abs(S[10000] - PROD) < 0.01, f"S_(h<=10^4) = {S[10000]:.5f} vs prod_(p<=10^7) = {PROD:.5f} (|diff| {abs(S[10000] - PROD):.1e} < 0.01)")
check(abs(PROD - math.pi / 4 * c_rep) < 1e-3, f"product to 10^7 vs (pi/4) c_exact: diff {abs(PROD - math.pi / 4 * c_rep):.1e}")

ZONES = [("axis", (1, 0)), ("3+4i", (3, 4)), ("4+3i", (4, 3)), ("5+12i", (5, 12)), ("12+5i", (12, 5)),
         ("8+15i", (8, 15)), ("15+8i", (15, 8)), ("7+24i", (7, 24)), ("24+7i", (24, 7)),
         ("20+21i", (20, 21)), ("21+20i", (21, 20))]

# all rational unit vectors with h' <= HN (for the overlap model)
HN = 1000
CUTS = (100, 300, 1000)
RN = math.isqrt(HN) + 1
allrays = {}
for a in range(-RN, RN + 1):
    for b in range(-RN, RN + 1):
        n = a * a + b * b
        if n == 0 or n > HN or n % 2 == 0 or gcd(a, b) != 1:
            continue
        w = complex(a, b) ** 2
        for e in (1, 1j, -1, -1j):
            ww = e * w
            allrays[(round(ww.real), round(ww.imag))] = n
for e in ((1, 0), (0, 1), (-1, 0), (0, -1)):
    allrays[e] = 1
Acoef = {}
for (wx, wy), h in allrays.items():
    mu = mobius(h)
    if mu:
        Acoef[(wx, wy)] = (4 / math.pi) * mu / phi(h)

reach = 2.4 * M / math.sqrt(math.sqrt(V_LO * X))   # max angular half-width of a zone (rad)
results = []
for name, (wx, wy) in ZONES:
    h = math.isqrt(wx * wx + wy * wy)
    phu = math.atan2(wy, wx)
    t = (A * wx + B * wy) / h; s = (B * wx - A * wy) / h
    zone = (t > 0)
    W = np.zeros(len(A))
    W[zone] = V(N[zone] / X) * F(s[zone] / (M * np.sqrt(t[zone])))
    k = W != 0
    wl = LAM[k] * W[k]
    psi = {mm: np.exp(1j * mm * ARG[k]) for mm in (0, 4)}
    Z = {mm: complex(np.sum(wl * psi[mm] * ee(-RAD[k]))) for mm in (0, 4)}
    ctl = [complex(np.sum(wl * psi[4] * ee(-c * RAD[k]))) for c in CONTROLS] + \
          [complex(np.sum(wl * psi[0] * ee(-c * RAD[k]))) for c in CONTROLS]
    sd = math.sqrt(float(np.sum(wl ** 2)) / 2)
    inert_part = complex(np.sum((wl * psi[4] * ee(-RAD[k]))[inert[k]]))
    # theorem main term and overlap model
    mu_h = mobius(h)
    Ah = (4 / math.pi) * mu_h / phi(h)
    u = complex(wx, wy) / h
    own = {mm: Ah * u ** mm * (I_q4 if mm == 4 else I_q0) for mm in (0, 4)}
    nb = {c: {0: 0j, 4: 0j} for c in CUTS}         # neighbouring rays (all rational directions != u in reach)
    nrays = 0
    for (vx, vy), Av in Acoef.items():
        if (vx, vy) == (wx, wy):
            continue
        phx = math.atan2(vy, vx)
        d = (phx - phu + math.pi) % TWO_PI - math.pi
        if abs(d) > reach + 0.01:
            continue
        nrays += 1
        Jv = Jm(phu, phx)
        for c in CUTS:
            if allrays[(vx, vy)] <= c:
                for mm in (0, 4):
                    nb[c][mm] += Av * Jv[mm]
    model = {mm: own[mm] + nb[HN][mm] for mm in (0, 4)}
    results.append((name, h, Ah, Z, own, model, nb, sd, ctl, inert_part, int(k.sum()), nrays))
    print(f"   zone {name:7s} h={h:3d} A_h={Ah:+.5f}: {int(k.sum())} primes, {nrays} rays (h'<={HN}) in reach ({time.time() - T0:.0f} s)")

print()
print("Z = zone sum over Gaussian prime elements without the inert primes q = 3 mod 4 (these lie on the axis; listed below)")
print("own = theorem's main term A_h u^m int G dA;  nb = heuristic contribution of the other rational rays in reach (h' <= 1000)")
print("zone      m      Z                        own (theorem)            nb (h'<=1000)        sd   (Z-own)/sd (Z-own-nb)/sd  nb cuts 100/300/1000 (in sd)")
for (name, h, Ah, Z, own, model, nb, sd, ctl, inert_part, npr, nrays) in results:
    for mm in (4, 0):
        Zm = Z[mm] - inert_part                     # inert primes: psi_m(q) = 1, e(-q) = 1, s = 0
        dth = abs(Zm - own[mm]) / sd
        dmo = abs(Zm - model[mm]) / sd
        cuts = " ".join(f"{abs(nb[c][mm]) / sd:4.2f}" for c in CUTS)
        print(f"{name:8s} {mm:2d}  {Zm.real:9.0f} {Zm.imag:+9.0f}i   {own[mm].real:9.0f} {own[mm].imag:+9.0f}i   "
              f"{nb[HN][mm].real:8.0f} {nb[HN][mm].imag:+8.0f}i  {sd:5.0f}  {dth:6.2f}  {dmo:6.2f}       {cuts}")
print()
print("inert primes (z = q = 3 mod 4, on the axis) in the axis zone: %.0f (real, same for m = 0, 4; subtracted above)"
      % results[0][9].real)
ax = results[0]
for mm in (4, 0):
    full = ax[3][mm]
    print(f"   axis m={mm}: with the inert primes Z = {full.real:.0f} {full.imag:+.0f}i, (Z - own)/sd = {abs(full - ax[4][mm]) / ax[7]:.2f}")
allctl = np.array([c for r in results for c in r[8]])
sds = np.array([r[7] for r in results for _ in r[8]])
rms_re = math.sqrt(float(np.mean((allctl.real / sds) ** 2)))
rms_im = math.sqrt(float(np.mean((allctl.imag / sds) ** 2)))
print(f"controls (4 frequencies x 2 m x {len(results)} zones): rms/sd real {rms_re:.2f}, imag {rms_im:.2f}")
check(0.5 < rms_re < 2.0 and 0.5 < rms_im < 2.0, "control sums are noise-sized (rms/sd within [0.5, 2])")
check(abs(results[0][3][4] - results[0][9] - results[0][4][4]) <= 4 * results[0][7] and
      abs(results[0][3][0] - results[0][9] - results[0][4][0]) <= 4 * results[0][7],
      "axis zone: Z within 4 sd of the theorem's main term (unconditional case h = 1)")
for (name, h, Ah, Z, own, model, nb, sd, ctl, inert_part, npr, nrays) in results:
    for mm in (4, 0):
        Zm = Z[mm] - inert_part
        check(abs(Zm - model[mm]) <= 4 * sd, f"zone {name} m={mm}: |Z - (own + nb)| = {abs(Zm - model[mm]) / sd:.2f} sd <= 4 (heuristic overlap model)")
# pooled amplitude of the theorem's own terms for the rays with A_h != 0 (m = 4), regression of Z on the model
num = den = 0j
for (name, h, Ah, Z, own, model, nb, sd, ctl, inert_part, npr, nrays) in results:
    if name == "axis" or Ah == 0:
        continue
    Zm = Z[4] - inert_part
    num += (Zm.conjugate() * model[4]).real / sd ** 2
    den += abs(model[4]) ** 2 / sd ** 2
amp = (num / den).real
amp_err = 1 / math.sqrt(den.real)
print(f"pooled amplitude of the Pythagorean zones (m=4, Z regressed on the overlap model): {amp:.4f} +- {amp_err:.4f}")
check(abs(amp - 1) < 4 * amp_err, "pooled Pythagorean amplitude within 4 errors of 1")

# ================================================================== F. approach to the theorem's regime (model only)
print("== F. neighbour cross-talk |nb|/|own| (h' <= 1000) as X grows, M = 3, same V (no primes involved)")
for name, (wx, wy) in (("axis", (1, 0)), ("3+4i", (3, 4)), ("5+12i", (5, 12)), ("20+21i", (20, 21)), ("7+24i", (7, 24))):
    h = math.isqrt(wx * wx + wy * wy); phu = math.atan2(wy, wx)
    Ah = (4 / math.pi) * mobius(h) / phi(h)
    row = []
    for Xv in (10 ** 7, 10 ** 9, 10 ** 11, 10 ** 13):
        rch = 2.4 * M / (V_LO * Xv) ** 0.25
        own4 = Ah * (complex(wx, wy) / h) ** 4 * own_main(4, Xv)
        nb4 = 0j; cnt = 0
        for (vx, vy), Av in Acoef.items():
            if (vx, vy) == (wx, wy):
                continue
            phx = math.atan2(vy, vx)
            d = (phx - phu + math.pi) % TWO_PI - math.pi
            if abs(d) > rch + 0.01:
                continue
            cnt += 1
            nb4 += Av * Jm(phu, phx, Xv)[4]
        row.append(f"X=1e{int(round(math.log10(Xv)))}: {cnt:3d} rays, |nb|/X^(3/4) = {abs(nb4) / Xv ** 0.75:.4f}"
                   + (f", |nb|/|own| = {abs(nb4) / abs(own4):.3f}" if Ah else ""))
    print(f"   {name:7s} " + " | ".join(row))

# ================================================================== E. elementary statement
print("== E. Re(z conj w) mod h for Gaussian primes in the flat part of the zone (|s| <= M sqrt t, N in [0.05X, X])")
for name, (wx, wy) in ZONES:
    h = math.isqrt(wx * wx + wy * wy)
    if h == 1:
        continue
    t = (A * wx + B * wy) / h; s = (B * wx - A * wy) / h
    flat = (t > 0) & (np.abs(s) <= M * np.sqrt(np.maximum(t, 0))) & (N >= V_LO * X)
    r = (A[flat] * wx + B[flat] * wy) % h
    lam = LAM[flat]
    bad = int(np.sum(np.gcd(r, h) != 1))
    cnt = np.bincount(r, minlength=h)
    units = [x for x in range(h) if gcd(x, h) == 1]
    obs = cnt[units]
    exp_ = obs.sum() / len(units)
    chi2 = float(np.sum((obs - exp_) ** 2 / exp_))
    dof = len(units) - 1
    zchi = (chi2 - dof) / math.sqrt(2 * dof)
    mean = complex(np.sum(lam * ee(-r / h)) / np.sum(lam))
    pred = mobius(h) / phi(h)
    sdm = math.sqrt(float(np.sum(lam ** 2)) / 2) / float(np.sum(lam))
    print(f"   {name:7s} h={h:3d}: {int(flat.sum())} primes, residues not coprime to h: {bad}; chi2 = {chi2:.1f} on {dof} dof"
          f" (z = {zchi:+.2f}); mean e(-r/h) = {mean.real:+.4f} {mean.imag:+.4f}i vs mu/phi = {pred:+.4f} (sd {sdm:.4f})")
    check(bad == 0, f"{name}: Re(z conj w) is never 0 mod a prime factor of h")
    check(abs(zchi) < 4, f"{name}: uniform on (Z/h)^* (chi2 z-score {zchi:+.2f})")
    check(abs(mean - pred) < 4 * sdm, f"{name}: mean e(-r/h) within 4 sd of mu(h)/phi(h)")
# lattice-point control for 3 + 4i (all lattice points in the same flat zone)
cnt = np.zeros(5, np.int64)
for a in range(1, math.isqrt(X) + 1):
    Bm = math.isqrt(X - a * a)
    b = np.arange(-Bm, Bm + 1, dtype=np.int64)
    n = a * a + b * b
    t = (3 * a + 4 * b) / 5; s = (3 * b - 4 * a) / 5
    sel = (t > 0) & (np.abs(s) <= M * np.sqrt(np.maximum(t, 0))) & (n >= V_LO * X)
    cnt += np.bincount((3 * a + 4 * b[sel]) % 5, minlength=5)
share = cnt / cnt.sum()
print(f"   lattice points in the 3+4i flat zone: residues 0..4 of 3a+4b mod 5: shares {np.round(share, 5)}")
check(np.max(np.abs(share - 0.2)) < 0.002, "lattice points: 3a + 4b mod 5 uniform on all five residues (including 0)")

print(f"time {time.time() - T0:.0f} s")
if FAILS:
    print(f"{len(FAILS)} FAIL(s)")
    raise AssertionError(FAILS[0])
print("ALL CHECKS PASSED")
