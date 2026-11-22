"""Round 89 (preregister/round89.txt): the k-th neighbour law for twin primes. Round 59's twin-gap law is the z^0
coefficient of G_a(z) = exp(sum_h ln(1 + p_h (z-1)) + 1/2 sum_{h!=k} R_hk w_h w_k) over the interior twin starts h of
class a; the z^{k-1} coefficient gives the distance to the k-th next twin pair. Data: twin pairs in [1e9, 2e9).
usage: r89.py"""
import os, json, math, time
import numpy as np
from scipy import stats
import common as C
os.environ.setdefault('R59_Y', '13')
import r59

OUT = r59.OUT
DMAX = 1200; NZ = 32; KS = (1, 2, 3)
Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)
L0, L1 = 10 ** 9, 2 * 10 ** 9
LAMS = [math.log(L0), math.log((L0 + L1) / 2), math.log(L1)]
BANDS = [(lo, lo + 24) for lo in range(6, 600, 30)] + [(lo, lo + 54) for lo in range(606, 1200, 60)]


def law():
    out = {k: np.zeros((len(LAMS), DMAX // 6)) for k in KS}
    lnS02 = sum(math.log(1 - 2 / q) - 2 * math.log(1 - 1 / q) for q in r59.MID) + r59.TAIL[2]
    for d in range(6, DMAX + 1, 6):
        h, lnB, lns, lnr = r59.lnS_parts(d)
        cls = r59.CLS[(np.gcd(r59.CLS + d, r59.PP) == 1) & (np.gcd(r59.CLS + d + 2, r59.PP) == 1)]
        if len(cls) == 0:
            continue
        M = ((np.gcd(cls[:, None] + h[None, :], r59.PP) == 1) & (np.gcd(cls[:, None] + h[None, :] + 2, r59.PP) == 1)).astype(float) if len(h) else np.zeros((len(cls), 0))
        R = np.expm1(lnr)
        for li, lam in enumerate(LAMS):
            rho0 = 1 / (lam * r59.DEN); p = rho0 ** 2 * np.exp(lns)
            logG = np.zeros((len(cls), NZ), complex)
            if len(h):
                logG += M @ np.log1p(p[:, None] * (Z[None, :] - 1))
            if len(h) > 1:
                W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))
                for m in range(NZ):
                    X = M * W[:, m][None, :]
                    logG[:, m] += 0.5 * np.sum((X.real @ R + 1j * (X.imag @ R)) * X, axis=1)
            coef = np.fft.fft(np.exp(logG), axis=1).real / NZ
            for k in KS:
                out[k][li, d // 6 - 1] = rho0 ** 2 * math.exp(lnB - lnS02) * coef[:, k - 1].sum() / len(r59.CLS)
    # average over the window as round 59 does (weights e^u/u^2 over ln x, linear interpolation between the 3 lam)
    u = np.linspace(LAMS[0], LAMS[-1], 401); wt = np.exp(u) / u ** 2
    avg = {}
    for k in KS:
        A = out[k]; avg[k] = np.array([np.sum(np.interp(u, LAMS, A[:, j]) * wt) / np.sum(wt) for j in range(A.shape[1])])
    return avg


def data():
    cnt = {k: {} for k in KS}; carry = np.zeros(0, np.int64); lastP = None
    for s0 in range(L0, L1, 10 ** 8):
        P = np.nonzero(C.primes_mask(s0, 10 ** 8))[0].astype(np.int64) + s0
        if lastP is not None:
            P = np.concatenate([[lastP], P])
        lastP = int(P[-1])
        tw = P[:-1][np.diff(P) == 2]
        T = np.concatenate([carry, tw])
        for k in KS:
            d = T[k:] - T[:-k]
            if len(carry):
                d = d[max(0, len(carry) - k):]
            v, c = np.unique(d[d <= DMAX], return_counts=True)
            for a, b in zip(v.tolist(), c.tolist()):
                cnt[k][a] = cnt[k].get(a, 0) + b
        carry = T[-max(KS):]
    return cnt


def chi2(cnt, lawd):
    N = sum(cnt.values()); lawd = lawd / lawd.sum(); ds = np.arange(6, DMAX + 1, 6)
    x = 0.0; dof = 0; zmax = 0.0
    for lo, hi in BANDS:
        c = sum(v for dd, v in cnt.items() if lo <= dd <= hi); e = N * lawd[(ds >= lo) & (ds <= hi)].sum()
        if e >= 50:
            x += (c - e) ** 2 / e; dof += 1; zmax = max(zmax, abs(c - e) / math.sqrt(e))
    return x, dof - 1, zmax


if __name__ == '__main__':
    t0 = time.time()
    avg = law(); print('law', round(time.time() - t0), 's', flush=True)
    cnt = data(); print('data', round(time.time() - t0), 's', {k: sum(v.values()) for k, v in cnt.items()}, flush=True)
    res = {'lams': LAMS, 'DMAX': DMAX, 'k': {}}
    one = avg[1] / avg[1].sum(); ds = np.arange(6, DMAX + 1, 6)
    for k in KS:
        x, dof, zmax = chi2(cnt[k], avg[k])
        full = one.copy()
        for _ in range(k - 1):
            full = np.convolve(full, one)
        dd = 6 * (np.arange(len(full)) + k); ind = np.array([full[dd == d].sum() for d in ds])
        xi, dofi, zi = chi2(cnt[k], ind)
        res['k'][str(k)] = {'pairs': int(sum(cnt[k].values())), 'law_mass_within_DMAX': float(avg[k].sum()),
                            'chi2': x, 'dof': dof, 'q999': float(stats.chi2.ppf(0.999, dof)), 'max_abs_z': zmax,
                            'independent_chi2': xi, 'independent_max_abs_z': zi}
        print(k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in res['k'][str(k)].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'r89.json'), 'w'), indent=1)
