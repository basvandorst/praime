#!/usr/bin/env python3
"""mock_test.py -- predictions and readings for the cumulants of the low-zero linear statistic over all chi mod q.

Statistic (from mock_sum.py): P_s(chi) = W_f(chi) - Wbar_a(chi), fhat = (1-|u|/s)_+, over the q-2 non-principal chi.
Cumulants kappa_m (m = 2, 3, 4) are population cumulants over the family (it is the whole family, not a sample).

Per (q, s), frozen predictions (mode 'predict', computed from primes only, no character values):
  St_m   independent phases (Steinhaus) cumulants, prime by prime  [the "diagonal": exact equality solutions]
  R_m    model-R ("no-wrap") deficit of the cumulants: tuples with 0 < |prod_S n - prod_S' n| < q cannot be congruent
         (weight -1/(q-2) each, exact), tuples with |..| >= q are "random" (expected weight 0)   [derivation.md s.4]
  U_m    random-matrix (U(N), N -> infinity) excess: m = 2: -int (|u|-1)_+ fhat^2; m = 3, 4: kappa_m^U(fhat)
  G_m    Gaussian: kappa_2 = int min(|u|,1) fhat^2, kappa_3 = kappa_4 = 0
Readings (mode 'read'): E_m = kappa_m - St_m (excess over the diagonal), S_m = E_m - R_m (the wrap-around,
"structured" part), compared with U_m; se from the delta method (moments to order 8, n_eff = (q-1)/2 conjugate pairs).

usage:
  mock_test.py predict --q 30011,30013 --sigmas 0.3,...,1.4 --out predictions.json
  mock_test.py read --pred predictions.json --data DIR       (DIR holds q<q>.npz from mock_sum.py)
  mock_test.py pilot --data pilot/q10007.npz                 (adds a bootstrap check of the delta-method se)
Self-contained: numpy only.
"""
import argparse, hashlib, json, math, os, sys, time
import numpy as np

# ---------------------------------------------------------------- primes
def simple_sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False; s[4::2] = False
    for i in range(3, math.isqrt(n) + 1, 2):
        if s[i]:
            s[i * i::2 * i] = False
    return np.nonzero(s)[0].astype(np.int64)


def primes_upto(X, seg=10_000_000):
    base = simple_sieve(math.isqrt(X) + 1)
    out = []
    lo = 0
    while lo <= X:
        hi = min(X + 1, lo + seg)
        s = np.ones(hi - lo, dtype=bool)
        if lo == 0:
            s[:min(2, hi)] = False
        for p in base:
            p = int(p)
            if p * p >= hi:
                break
            st = max(p * p, ((lo + p - 1) // p) * p)
            s[st - lo::p] = False
        out.append(lo + np.nonzero(s)[0].astype(np.int64))
        lo = hi
    return np.concatenate(out), base


def prime_power_list(q, X):
    P, base = primes_upto(X)
    P = P[P != q]
    ns, ws, isp = [P], [np.log(P.astype(float))], [np.ones(len(P), dtype=bool)]
    for p in base:
        p = int(p)
        if p == q or p * p > X:
            continue
        k = 2
        while p ** k <= X:
            ns.append(np.array([p ** k], dtype=np.int64)); ws.append(np.array([math.log(p)]))
            isp.append(np.array([False])); k += 1
    n = np.concatenate(ns); lam = np.concatenate(ws); ip = np.concatenate(isp)
    o = np.argsort(n)
    return n[o], lam[o], ip[o]


def weights(nn, lam, q, s):
    L = math.log(q)
    return lam / np.sqrt(nn) * np.clip(1.0 - np.log(nn) / (L * s), 0.0, None)

# ---------------------------------------------------------------- Steinhaus (diagonal)
def one_prime(coef):
    K = len(coef); M = 4 * K + 4
    th = 2 * np.pi * np.arange(M) / M
    X = sum(c * np.cos((k + 1) * th) for k, c in enumerate(coef))
    m2 = np.mean(X ** 2)
    return m2, np.mean(X ** 3), np.mean(X ** 4) - 3 * m2 * m2


def steinhaus(nn, lam, ip, q, s):
    L = math.log(q); X = q ** s
    keep = nn < X
    nn, lam, ip = nn[keep], lam[keep], ip[keep]
    a = weights(nn, lam, q, s)
    root = math.isqrt(int(X))
    big = ip & (nn > root)
    A = (2 / L) * a[big]
    k2 = np.sum(A * A) / 2; k3 = 0.0; k4 = -3 / 8 * np.sum(A ** 4)
    for p in nn[ip & (nn <= root)]:
        p = int(p); coef = []; k = 1
        while p ** k < X:
            coef.append(-(2 / L) * math.log(p) / p ** (k / 2) * (1 - k * math.log(p) / (L * s))); k += 1
        c2, c3, c4 = one_prime(coef)
        k2 += c2; k3 += c3; k4 += c4
    return float(k2), float(k3), float(k4)

# ---------------------------------------------------------------- model R (no-wrap deficit)
class _Win:
    def __init__(self, nn, aa, q):
        self.nn, self.aa, self.q = nn, aa, q
        self.C = np.concatenate([[0.0], np.cumsum(aa)])

    def le(self, x):
        return self.C[np.searchsorted(self.nn, x, side='right')]

    def near(self, x):
        nn, aa, C, q = self.nn, self.aa, self.C, self.q
        x = np.asarray(x, dtype=np.int64)
        idx = np.searchsorted(nn, x)
        eq = np.zeros(len(x)); ok = idx < len(nn)
        hit = np.zeros(len(x), dtype=bool); hit[ok] = nn[idx[ok]] == x[ok]
        eq[hit] = aa[idx[hit]]
        return C[np.searchsorted(nn, x + q, side='left')] - C[np.searchsorted(nn, x - q, side='right')] - eq


def _pairs(nn, aa, Z):
    ps, ws = [], []
    for i in range(len(nn)):
        n1 = int(nn[i])
        if n1 * n1 > Z:
            break
        j = np.searchsorted(nn, Z // n1, side='right')
        if j <= i:
            continue
        ps.append(n1 * nn[i:j]); w = aa[i] * aa[i:j] * 2.0; w[0] = aa[i] * aa[i]; ws.append(w)
    if not ps:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    p = np.concatenate(ps); w = np.concatenate(ws); o = np.argsort(p)
    return p[o], w[o]


def model_R(nn_all, lam_all, q, s, mmax=4, ycut_factor=50.0, delta=2e-3):
    L = math.log(q)
    keep = nn_all < q ** s
    nn = nn_all[keep]; aa = weights(nn, lam_all[keep], q, s)
    W = _Win(nn, aa, q); Xs = int(nn[-1]); lim = Xs + q
    N = {'10': float(W.le(q))}
    small = nn[nn <= q]
    N['20'] = float(sum(aa[i] * W.le(q // int(x)) for i, x in enumerate(small)))
    N['11'] = float(np.sum(aa * W.near(nn)))
    p2, w2 = _pairs(nn, aa, q)
    N['30'] = float(sum(w * W.le(q // int(x)) for x, w in zip(p2, w2)))
    t21 = 0.0
    for i in range(len(nn)):
        n1 = int(nn[i])
        if n1 * n1 >= lim:
            break
        j = np.searchsorted(nn, (lim - 1) // n1, side='right')
        if j <= i:
            continue
        w = aa[i] * aa[i:j] * 2.0; w[0] = aa[i] * aa[i]
        t21 += float(np.sum(w * W.near(n1 * nn[i:j])))
    N['21'] = t21
    if mmax >= 4:
        N['40'] = float(sum(w * np.sum(w2[p2 <= q // int(x)]) for x, w in zip(p2, w2)))
        pz, wz = _pairs(nn, aa, (lim - 1) // 2)
        t31 = 0.0
        for i3 in range(len(nn)):
            n3 = int(nn[i3]); j = np.searchsorted(pz, (lim - 1) // n3, side='right')
            if j == 0:
                break
            t31 += aa[i3] * float(np.sum(wz[:j] * W.near(pz[:j] * n3)))
        N['31'] = t31
        del pz, wz
        Ycut = int(ycut_factor * q)
        pe, we = _pairs(nn, aa, Ycut + 2 * q)
        Ce = np.concatenate([[0.0], np.cumsum(we)])
        _, start, cnt = np.unique(pe, return_index=True, return_counts=True)
        eq_each = np.repeat(np.add.reduceat(we, start), cnt) if len(we) else np.zeros(0)
        A = pe < Ycut
        lo = np.searchsorted(pe, pe[A] - q, side='right'); hi = np.searchsorted(pe, pe[A] + q, side='left')
        exact = float(np.sum(we[A] * (Ce[hi] - Ce[lo] - eq_each[A])))
        del pe, we, Ce, eq_each, lo, hi
        nb = int((2 * math.log(Xs) + 2 * delta) / delta) + 2
        H1 = np.bincount((np.log(nn) / delta).astype(np.int64), weights=aa, minlength=nb)[:nb]
        H2 = np.convolve(H1, H1)
        xk = np.exp((np.arange(len(H2)) + 1.0) * delta)
        sel = xk >= Ycut
        smooth = float(2 * q * np.sum(H2[sel] ** 2 / (xk[sel] * delta)))
        N['22'] = exact + smooth; N['22_smooth_share'] = smooth / max(exact + smooth, 1e-300)
    c = -1.0 / (q - 2)
    Dm = {1: (-1 / L) * c * 2 * N['10'], 2: (1 / L ** 2) * c * (2 * N['20'] + 2 * N['11']),
          3: (-1 / L ** 3) * c * (2 * N['30'] + 6 * N['21'])}
    if mmax >= 4:
        Dm[4] = (1 / L ** 4) * c * (2 * N['40'] + 8 * N['31'] + 6 * N['22'])
    return Dm, N


def R_cumulants(Dm, st):
    k2s, k3s, k4s = st
    mu1 = Dm[1]; mu2 = k2s + Dm[2]; mu3 = k3s + Dm[3]
    out = {'R1': mu1, 'R2': mu2 - mu1 ** 2 - k2s, 'R3': mu3 - 3 * mu2 * mu1 + 2 * mu1 ** 3 - k3s}
    if 4 in Dm:
        mu4 = k4s + 3 * k2s ** 2 + Dm[4]
        out['R4'] = mu4 - 4 * mu3 * mu1 - 3 * mu2 ** 2 + 12 * mu2 * mu1 ** 2 - 6 * mu1 ** 4 - k4s
    return out

# ---------------------------------------------------------------- random matrix theory
def fhat_tri(u, s):
    return np.clip(1.0 - np.abs(u) / s, 0.0, None)


def toeplitz_cumulants(N, s):
    """exact cumulants 2..4 of Z = sum_n (1/N) fhat(n/N) Tr U^n over U(N) (Heine). N = 300 equals N = inf to 1e-5."""
    K = int(np.floor(s * N)); n = np.arange(-K, K + 1)
    c = fhat_tri(n / N, s) / N
    pw = [None, c]
    for _ in range(3):
        pw.append(np.convolve(pw[-1], c))
    idx = np.arange(N)[:, None] - np.arange(N)[None, :]
    def T(co):
        Lh = (len(co) - 1) // 2; M = np.zeros((N, N)); ok = np.abs(idx) <= Lh; M[ok] = co[idx[ok] + Lh]; return M
    A1, A2, A3, A4 = (T(pw[k]) for k in range(1, 5))
    A11 = A1 @ A1; tr = np.trace
    k2 = tr(A2) - tr(A11)
    k3 = tr(A3) - 3 * tr(A1 @ A2) + 2 * tr(A11 @ A1)
    k4 = tr(A4) - 4 * tr(A1 @ A3) - 3 * tr(A2 @ A2) + 12 * tr(A11 @ A2) - 6 * tr(A11 @ A11)
    return float(k2), float(k3), float(k4)


def rmt_row(s, N=300):
    k2, k3, k4 = toeplitz_cumulants(N, s)
    diag_inf = 2 * s * s / 12.0                                    # int |u| fhat^2 = s^2/6
    E2 = -2 * ((s - 1) ** 4 / 12 / s ** 2) if s > 1 else 0.0       # -int (|u|-1)_+ fhat^2, exact for the triangle
    return {'U2': E2, 'U3': k3, 'U4': k4, 'k2U': k2, 'G2': min(diag_inf, diag_inf + E2), 'diag_inf': diag_inf}

# ---------------------------------------------------------------- measured cumulants and se
def central_moments(ps):
    """ps[e] = sum x^e, e = 0..8 -> mean and central moments c_0..c_8."""
    n = ps[0]; m = ps[1] / n
    raw = ps / n
    c = np.zeros(9)
    for r in range(9):
        c[r] = sum(math.comb(r, j) * raw[j] * (-m) ** (r - j) for j in range(r + 1))
    return m, c


def cumulants_se(ps, n_eff):
    m, c = central_moments(ps)
    k2 = c[2]; k3 = c[3]; k4 = c[4] - 3 * c[2] ** 2
    S = np.array([[c[r + t] - c[r] * c[t] for t in range(1, 5)] for r in range(1, 5)])
    g = {2: np.array([0, 1, 0, 0.0]), 3: np.array([-3 * c[2], 0, 1, 0.0]), 4: np.array([-4 * c[3], -6 * c[2], 0, 1.0])}
    se = {k: math.sqrt(max(g[k] @ S @ g[k], 0.0) / n_eff) for k in g}
    return m, {2: k2, 3: k3, 4: k4}, se

# ---------------------------------------------------------------- modes
def do_predict(a):
    qs = [int(x) for x in a.q.split(',')]
    sig = [float(x) for x in a.sigmas.split(',')]
    res = json.load(open(a.out)) if os.path.exists(a.out) else {'cells': {}}
    rmt = {}
    for s in sig:
        rmt[s] = rmt_row(s)
    for q in qs:
        smax = max(sig)
        nn, lam, ip = prime_power_list(q, int(q ** smax) + 1)
        for s in sig:
            key = f'{q}_{s:.2f}'
            if key in res['cells']:
                continue
            t0 = time.time()
            st = steinhaus(nn, lam, ip, q, s)
            m4 = 4 if s <= a.s4max + 1e-9 else 3
            Dm, N = model_R(nn, lam, q, s, mmax=m4, ycut_factor=a.ycut)
            R = R_cumulants(Dm, st)
            cell = {'q': q, 's': s, 'logq': math.log(q), 'St2': st[0], 'St3': st[1], 'St4': st[2], **R, **rmt[s],
                    'N22_smooth_share': N.get('22_smooth_share', None)}
            res['cells'][key] = cell
            print(f'{key}: St {st[0]:.5f} {st[1]:+.6f} {st[2]:+.6f} | R {R["R2"]:+.6f} {R["R3"]:+.6f} '
                  f'{R.get("R4", float("nan")):+.6f} | U {rmt[s]["U2"]:+.6f} {rmt[s]["U3"]:+.6f} {rmt[s]["U4"]:+.6f}'
                  f' ({time.time()-t0:.0f}s)', flush=True)
            json.dump(res, open(a.out, 'w'), indent=1, sort_keys=True)
        del nn, lam, ip
    blob = json.dumps(res['cells'], sort_keys=True).encode()
    print('sha256(cells) =', hashlib.sha256(blob).hexdigest())


def read_one(path, pred=None):
    z = np.load(path)
    q = int(z['q']); sig = z['sig']; n_eff = (q - 1) / 2
    rows = []
    for j, s in enumerate(sig):
        m, k, se = cumulants_se(z['psums'][j, 0], n_eff)
        st = z['diag'][j]
        row = {'q': q, 's': float(s), 'mean': m, 'k2': k[2], 'k3': k[3], 'k4': k[4], 'se2': se[2], 'se3': se[3],
               'se4': se[4], 'St2': st[0], 'St3': st[1], 'St4': st[2]}
        for mm in (2, 3, 4):
            row[f'E{mm}'] = k[mm] - st[mm - 2]
        if pred is not None:
            c = pred['cells'].get(f'{q}_{float(s):.2f}')
            if c is not None:
                row.update({kk: c[kk] for kk in ('R1', 'R2', 'R3', 'R4', 'U2', 'U3', 'U4') if kk in c})
                assert abs(c['St2'] - st[0]) < 1e-9 * max(1, st[0]) + 1e-12, (c['St2'], st[0])
                for mm in (2, 3, 4):
                    if f'R{mm}' in row:
                        row[f'S{mm}'] = row[f'E{mm}'] - row[f'R{mm}']
        rows.append(row)
    return rows


def do_read(a):
    pred = json.load(open(a.pred))
    files = sorted(f for f in os.listdir(a.data) if f.startswith('q') and f.endswith('.npz'))
    allrows = []
    for f in files:
        allrows += read_one(os.path.join(a.data, f), pred)
    hdr = ('q        s    k2       St2      E2        R2        S2       U2       se2    | k3        St3       E3'
           '        R3        S3       U3       se3    | k4        St4       E4        R4        S4       U4       se4')
    print(hdr)
    for r in allrows:
        def g(k):
            return r.get(k, float('nan'))
        print(f"{r['q']:<8d} {r['s']:.2f} {r['k2']:.5f} {r['St2']:.5f} {r['E2']:+.6f} {g('R2'):+.6f} {g('S2'):+.6f} "
              f"{g('U2'):+.6f} {r['se2']:.5f} | {r['k3']:+.6f} {r['St3']:+.6f} {r['E3']:+.6f} {g('R3'):+.6f} "
              f"{g('S3'):+.6f} {g('U3'):+.6f} {r['se3']:.5f} | {r['k4']:+.6f} {r['St4']:+.6f} {r['E4']:+.6f} "
              f"{g('R4'):+.6f} {g('S4'):+.6f} {g('U4'):+.6f} {r['se4']:.5f}")
    # readings
    print('\nREADINGS (z = (S - U)/se per cell; pooled over moduli with inverse-variance weights)')
    for name, m, sl in (('T0 m=2', 2, (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), ('T0 m=3', 3, (0.3, 0.4, 0.5, 0.6)),
                        ('T0 m=4', 4, (0.3, 0.4)), ('T1 m=2', 2, (1.0,)), ('T1 m=3', 3, (0.7, 0.8, 0.9, 1.0)),
                        ('T1 m=4', 4, (0.5, 0.6, 0.7)),
                        ('T2 m=2', 2, (1.1, 1.2, 1.3, 1.4)), ('T2 m=3', 3, (1.1, 1.2, 1.3, 1.4)),
                        ('T2 m=4', 4, (0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4))):
        for s in sl:
            rr = [r for r in allrows if abs(r['s'] - s) < 1e-9 and f'S{m}' in r]
            if not rr:
                continue
            S = np.array([r[f'S{m}'] for r in rr]); U = np.array([r[f'U{m}'] for r in rr])
            se = np.array([r[f'se{m}'] for r in rr]); w = 1 / se ** 2
            z = (S - U) / se
            Sp = np.sum(w * S) / np.sum(w); sep = 1 / math.sqrt(np.sum(w))
            chi2 = float(np.sum(z ** 2))
            ratio = Sp / U.mean() if abs(U.mean()) > 0 else float('nan')
            print(f'{name} s={s:.1f}: n_q={len(rr)} pooled S={Sp:+.6f} +- {sep:.6f}  U={U.mean():+.6f}  '
                  f'z_pool={(Sp-U.mean())/sep:+.2f}  z0={Sp/sep:+.2f}  chi2(S=U)={chi2:.1f}/{len(rr)}  S/U={ratio:+.3f}')
            if name.startswith('T1'):
                for gname, lo, hi in (('F3', 2e4, 5e4), ('F5', 5e4, 2e5), ('F35', 2e5, 5e5), ('F6', 5e5, 5e6)):
                    g = [r for r in rr if lo < r['q'] < hi]
                    if g:
                        wg = np.array([1 / r[f'se{m}'] ** 2 for r in g])
                        Sg = np.sum(wg * np.array([r[f'S{m}'] for r in g])) / wg.sum()
                        print(f'      group {gname}: S={Sg:+.6f} +- {1/math.sqrt(wg.sum()):.6f}  S*log q={Sg*np.mean([math.log(r["q"]) for r in g]):+.5f}')
    # noise check: within-group scatter of S over the T0 cells
    print('\nNOISE CHECK (T0 cells): chi2 of (S - group mean)/se within groups')
    for gname, lo, hi in (('F3', 2e4, 5e4), ('F5', 5e4, 2e5), ('F35', 2e5, 5e5), ('F6', 5e5, 5e6)):
        chi2 = 0.0; dof = 0
        for m, sl in ((2, (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), (3, (0.3, 0.4, 0.5, 0.6)), (4, (0.3, 0.4))):
            for s in sl:
                g = [r for r in allrows if lo < r['q'] < hi and abs(r['s'] - s) < 1e-9 and f'S{m}' in r]
                if len(g) < 2:
                    continue
                x = np.array([r[f'S{m}'] for r in g]); se = np.array([r[f'se{m}'] for r in g])
                chi2 += float(np.sum(((x - x.mean()) / se) ** 2)); dof += len(g) - 1
        if dof:
            print(f'  {gname}: chi2/dof = {chi2:.1f}/{dof}')
    # T3: extrapolated wrap-around fourth cumulant
    A = {1.2: 0.0037, 1.3: 0.0072, 1.4: 0.0159}
    print('\nT3 (m=4): weighted fit S_4 = alpha + beta/log q over moduli < 5e5')
    for s in (1.2, 1.3, 1.4):
        rr = [r for r in allrows if abs(r['s'] - s) < 1e-9 and 'S4' in r and r['q'] < 5e5]
        if len({r['q'] for r in rr}) < 3:
            print(f'  s={s}: fewer than 3 moduli, no fit')
            continue
        x = np.array([1 / math.log(r['q']) for r in rr]); y = np.array([r['S4'] for r in rr])
        w = np.array([1 / r['se4'] ** 2 for r in rr])
        M = np.vstack([np.ones_like(x), x]).T
        cov = np.linalg.inv(M.T @ (M * w[:, None])); coef = cov @ (M.T @ (w * y))
        sea = math.sqrt(cov[0, 0]); U4 = rr[0]['U4']
        c = pred['cells'][f"{rr[0]['q']}_{s:.2f}"]; Dp = c['Dpilot4']
        den = math.sqrt(sea ** 2 + A[s] ** 2)
        tc = (coef[0] - U4) / den; ts = (coef[0] - U4 - Dp) / den
        print(f'  s={s}: alpha={coef[0]:+.5f} +- {sea:.5f} (allowance {A[s]})  beta={coef[1]:+.4f}  U4={U4:+.5f}  '
              f't_conv={tc:+.2f} [{"pass" if abs(tc) <= 3 else "fail"}]  t_stall={ts:+.2f} '
              f'[{"rejected" if abs(ts) >= 4 else "not rejected"}]')


def do_pilot(a):
    z = np.load(a.data)
    q = int(z['q']); P = z['P']; sig = z['sig']
    rng = np.random.default_rng(1)
    k = np.arange(1, q - 1)
    pair_rep = k[(k < (q - 1) / 2)]            # one of each conjugate pair; plus the real character k=(q-1)/2
    leg = (q - 1) // 2
    print('s    k2       se2(delta) se2(boot) | k3        se3(delta) se3(boot) | k4        se4(delta) se4(boot)')
    for j, s in enumerate(sig):
        x = P[j]
        m, kk, se = cumulants_se(z['psums'][j, 0], (q - 1) / 2)
        xp = x[pair_rep - 1]; xl = x[leg - 1]
        bs = []
        for _ in range(300):
            b = rng.choice(xp, size=len(xp), replace=True)
            y = np.concatenate([b, b, [xl]])
            c = y - y.mean(); m2 = np.mean(c ** 2)
            bs.append((m2, np.mean(c ** 3), np.mean(c ** 4) - 3 * m2 ** 2))
        bs = np.array(bs).std(axis=0)
        print(f'{s:.1f}  {kk[2]:.5f}  {se[2]:.6f}  {bs[0]:.6f} | {kk[3]:+.6f}  {se[3]:.6f}  {bs[1]:.6f} | '
              f'{kk[4]:+.6f}  {se[4]:.6f}  {bs[2]:.6f}')


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='mode', required=True)
    p1 = sp.add_parser('predict'); p1.add_argument('--q', required=True); p1.add_argument('--sigmas', required=True)
    p1.add_argument('--out', required=True); p1.add_argument('--s4max', type=float, default=1.4)
    p1.add_argument('--ycut', type=float, default=50.0)
    p2 = sp.add_parser('read'); p2.add_argument('--pred', required=True); p2.add_argument('--data', required=True)
    p3 = sp.add_parser('pilot'); p3.add_argument('--data', required=True)
    a = ap.parse_args()
    {'predict': do_predict, 'read': do_read, 'pilot': do_pilot}[a.mode](a)


if __name__ == '__main__':
    main()
