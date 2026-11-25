"""r2_lattice.py (review, own code): what lattice geometry alone gives.
Part A (Q(i), m = 4) on N in (3e7, X1]: sum over lattice points z = a + bi, a >= 1, b >= 0 (one point per ideal, axis included)
or b >= 1 (open quadrant, as the prime sum), of e^{4 i arg z} e(-|z|) N^{-1/4} W(|z|) with weight kappa_P for points with
gcd(N, prod_{p<=P} p) = 1 (P = 0: all points; 2: a+b odd; 10, 100, 1000); kappa_P = (4/pi)/delta_P makes the mean density
equal to that of the ln-weighted Gaussian primes.  c_eff = sum / MT_1 (MT_1 = main term with c = 1).
Also the axis zone sums over the lattice (open quadrant, sharp N range), m = 0, 4, for the zone model.
Part C (C23) on N in (3e7, X1c]: lattice points of x^2+xy+6y^2 and 2x^2+xy+3y^2 with weights +1/2, -1/2, sieved the same way.
usage: r2_lattice.py A X1 | C X1"""
import os, sys, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/')
import r1_sums as R
SD = R.SD; TP = 2 * math.pi
PS = [10, 100, 1000]
SMALLP = [int(p) for p in R.small_primes(1000)]

def rough_flags(L, U):
    """for numbers L+1 .. U-1: dict P -> bool array (gcd(n, prod_{p<=P} p) == 1), P in PS (nested)."""
    n = U - L - 1; s = np.ones(n, bool); out = {}
    for P in PS:
        for p in SMALLP:
            if p > P: break
            if P > 10 and p <= (PS[PS.index(P) - 1]): continue
            st = (-(L + 1)) % p; s[st::p] = False
        out[P] = s.copy()
    return out

def hann(v, v0, v1):
    u = (v - v0) / (v1 - v0); return np.where((u > 0) & (u < 1), np.sin(np.pi * u) ** 2, 0.0)

def mt1(K1, kappa, Q, v0, v1, n=200001):
    """(3/2) K1 int W(v) (1 - i kappa / t) dv, t = 2 pi v / sqrt Q (degree 2, v = sqrt x)."""
    v = np.linspace(v0, v1, n); f = hann(v, v0, v1) * (1 - 1j * kappa / (TP * v / math.sqrt(Q)))
    h = v[1] - v[0]; return 1.5 * K1 * h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())

def grid_chunks(x0, x1, step=20_000_000):
    out = []; x = x0
    while x < x1: y = min(x1, x + step); out.append((x, y)); x = y
    return out

def partA(X1):
    v0, v1 = math.sqrt(3e7), math.sqrt(X1)
    acc = {k: 0j for k in ['all_ideal', 'all_open', 'odd_ideal', 'odd_open'] + ['P%d_ideal' % P for P in PS] + ['P%d_open' % P for P in PS]}
    zon = {'ax0_m0': 0j, 'ax0_m4': 0j, 'ax90_m4': 0j, 'n': 0}
    t0 = time.time()
    for (L, U) in grid_chunks(30_000_000, X1):
        rf = rough_flags(L, U)
        amax = int(R.isq(U - 1))
        for A0 in range(1, amax + 1, 300):
            A = np.arange(A0, min(A0 + 300, amax + 1), dtype=np.int64)
            lo = np.where(L + 1 - A * A > 0, R.isq(np.maximum(L - A * A, 0)) + 1, 0)
            hi = R.isq(np.maximum(U - 1 - A * A, -1)); hi = np.where(U - 1 - A * A >= 0, hi, -1)
            a, b = R.lattice_block(A, lo, hi)
            if len(a) == 0: continue
            N = a * a + b * b; m = (N > L) & (N < U); a, b, N = a[m], b[m], N[m]
            Nf = N.astype(np.float64); v = np.sqrt(Nf)
            z = (a + 1j * b.astype(np.float64)) / v; e4 = z ** 4; E = np.exp(-1j * TP * v)
            t = e4 * E * Nf ** -0.25 * hann(v, v0, v1)
            opn = b >= 1; odd = ((a + b) & 1) == 1
            acc['all_ideal'] += t.sum(); acc['all_open'] += t[opn].sum()
            acc['odd_ideal'] += t[odd].sum(); acc['odd_open'] += t[odd & opn].sum()
            idx = N - (L + 1)
            for P in PS:
                r = rf[P][idx]; acc['P%d_ideal' % P] += t[r].sum(); acc['P%d_open' % P] += t[r & opn].sum()
            # zones (sharp N range), octant a > b >= 1 and mirrors, all lattice points, weight 1
            mo = (a > b) & (b >= 1); ao = a[mo].astype(np.float64); bo = b[mo].astype(np.float64)
            u = bo / np.sqrt(ao); w = np.exp(-np.pi * u * u / 4); w[u > 12] = 0
            Eo = E[mo]; e4o = e4[mo]
            zon['ax0_m0'] += np.sum(w * Eo); zon['ax0_m4'] += np.sum(w * Eo * e4o); zon['ax90_m4'] += np.sum(w * Eo * np.conj(e4o))
        print('A (%d,%d] %.0fs' % (L, U, time.time() - t0), flush=True)
    return acc, zon, (v0, v1)

def form_points(Af, Bf, Cf, L, U):
    """all (x, y) with y >= 1 (any x) or y = 0, x >= 1, and L < f(x,y) < U."""
    D = Bf * Bf - 4 * Af * Cf
    ymax = math.isqrt(4 * Af * U // (-D)) + 2
    y = np.arange(0, ymax + 1, dtype=np.int64)
    dU = 4 * Af * (U - 1) + D * y * y; dL = 4 * Af * (L + 1) + D * y * y
    okU = dU >= 0; y = y[okU]; dU = dU[okU]; dL = dL[okU]
    sU = R.isq(dU) + 1; sL = np.where(dL > 0, R.isq(np.maximum(dL, 0)) - 1, -1)
    xs, ys = [], []
    for sgn in (1, -1):
        if sgn == 1: lo = (-Bf * y + np.maximum(sL, 0)) // (2 * Af) - 1; hi = (-Bf * y + sU) // (2 * Af) + 1
        else: lo = (-Bf * y - sU) // (2 * Af) - 1; hi = (-Bf * y - np.maximum(sL, 0)) // (2 * Af) + 1
        x, yy = R.lattice_block(y, lo, hi)      # lattice_block(A=y, blo, bhi) returns (y repeated, x values)
        xs.append(yy); ys.append(x)
    x = np.concatenate(xs); y = np.concatenate(ys)
    f = Af * x * x + Bf * x * y + Cf * y * y
    keep = (f > L) & (f < U) & ((y >= 1) | ((y == 0) & (x >= 1)))
    x, y, f = x[keep], y[keep], f[keep]
    # the two x-intervals may overlap when sL < 0: remove duplicates
    key = x * 1_000_003 + y; _, ui = np.unique(key, return_index=True)
    return x[ui], y[ui], f[ui]

def partC(X1):
    Q = 23; v0, v1 = math.sqrt(3e7), math.sqrt(X1)
    acc = {k: 0j for k in ['all'] + ['P%d' % P for P in PS]}
    t0 = time.time()
    for (L, U) in grid_chunks(30_000_000, X1, 5_000_000):
        rf = rough_flags(L, U)
        for (Af, Bf, Cf, sg) in ((1, 1, 6, 0.5), (2, 1, 3, -0.5)):
            x, y, f = form_points(Af, Bf, Cf, L, U)
            ff = f.astype(np.float64); v = np.sqrt(ff)
            t = 2 * sg * np.exp(-1j * TP * 2 * np.sqrt(ff / Q)) * ff ** -0.25 * hann(v, v0, v1)   # factor 2: (x,y) and (-x,-y)
            acc['all'] += t.sum(); idx = f - (L + 1)
            for P in PS: acc['P%d' % P] += t[rf[P][idx]].sum()
        print('C (%d,%d] %.0fs' % (L, U, time.time() - t0), flush=True)
    return acc, (v0, v1)

def delta_A(P):
    d = 1.0
    for p in SMALLP:
        if p > P: break
        if p == 2: d *= 0.5
        elif p % 4 == 1: d *= (1 - 1 / p) ** 2
        else: d *= 1 - 1 / p ** 2
    return d
def delta_C(P, D=23):
    d = 1.0
    for p in SMALLP:
        if p > P: break
        if p == D: d *= 1 - 1 / p
        else:
            ch = 1 if p == 2 else (1 if pow(p % D, (D - 1) // 2, D) == 1 else -1)
            rho = 2 * p - 1 if ch == 1 else 1
            d *= 1 - rho / p ** 2
    return d

if __name__ == '__main__':
    part = sys.argv[1]; X1 = int(float(sys.argv[2]))
    if part == 'A':
        acc, zon, (v0, v1) = partA(X1)
        K1 = (2 / 3) * complex(math.cos(-math.pi / 4), math.sin(-math.pi / 4)); M = mt1(K1, 1 / 16 - 4, 4, v0, v1)
        res = {'X1': X1, 'MT1': [M.real, M.imag], 'sums': {k: [v.real, v.imag] for k, v in acc.items()},
               'zones': {k: ([v.real, v.imag] if isinstance(v, complex) else v) for k, v in zon.items()}}
        kap = {'all': 4 / math.pi, 'odd': 8 / math.pi, **{'P%d' % P: (4 / math.pi) / delta_A(P) for P in PS}}
        res['c_eff'] = {}
        for k, v in acc.items():
            ce = kap[k.split('_')[0]] * v / M; res['c_eff'][k] = [ce.real, ce.imag]
            print('%-12s c_eff = %.5f %+.5fi' % (k, ce.real, ce.imag))
        json.dump(res, open(SD + '/r2_A.json', 'w'), indent=1)
    else:
        acc, (v0, v1) = partC(X1)
        K1 = (2 * math.sqrt(2) / 3) * 23 ** -0.25 * complex(math.cos(-math.pi / 4), math.sin(-math.pi / 4)); M = mt1(K1, 1 / 16, 23, v0, v1)
        L1 = 3 * math.pi / math.sqrt(23)
        res = {'X1': X1, 'MT1': [M.real, M.imag], 'c_eff': {}}
        for k, v in acc.items():
            kap = (1 / L1) / (1.0 if k == 'all' else delta_C(int(k[1:])))
            ce = kap * v / M; res['c_eff'][k] = [ce.real, ce.imag]; print('%-6s c_eff = %.5f %+.5fi' % (k, ce.real, ce.imag))
        json.dump(res, open(SD + '/r2_C.json', 'w'), indent=1)
