"""Round 83 (preregister/round83.txt): the distribution of the distance to the k-th neighbour, D_k = p_{n+k} - p_n,
from round 81's generating-function law (copied here), against the primes in [1e9, 1.1e9), for k = 2, 3, 4; compared
with Gallagher's Gamma(k) law and with the sum of k independent class-averaged gaps. usage: r83.py"""
import os, json, math
import numpy as np
from scipy import stats
import common as C
import r47

assert r47.Y == 13
OUT = r47.OUT
PP = r47.PP
LAM = math.log(1.05e9); RHO0 = 1 / (LAM * r47.RHO_DEN)
DMAX = 240; NZ = 48; KS = (1, 2, 3, 4)
Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)
NCLS = sum(1 for a in range(PP) if math.gcd(a, PP) == 1)


def law():
    dist = {k: np.zeros(DMAX // 2) for k in KS}
    for D in range(2, DMAX + 1, 2):
        st = r47.setup(D); M = st['M']; n = M.shape[1]; p = RHO0 * st['s1']
        logG = np.zeros((M.shape[0], NZ), complex)
        if n:
            logG += M @ np.log1p(p[:, None] * (Z[None, :] - 1))
        if n > 1:
            W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))
            for m in range(NZ):
                X = M * W[:, m][None, :]
                logG[:, m] += 0.5 * np.sum((X.real @ st['R'] + 1j * (X.imag @ st['R'])) * X, axis=1)
        coef = np.fft.fft(np.exp(logG), axis=1).real / NZ
        for k in KS:
            dist[k][D // 2 - 1] = RHO0 * st['S2'] * coef[:, k - 1].sum() / NCLS
    return dist


def chi2(counts, expect):
    keep = []; c = 0.0; e = 0.0; rows = []
    for ci, ei in zip(counts, expect):
        c += ci; e += ei
        if e >= 50:
            rows.append((c, e)); c = 0.0; e = 0.0
    if e > 0 and rows:
        rows[-1] = (rows[-1][0] + c, rows[-1][1] + e)
    x = sum((a - b) ** 2 / b for a, b in rows)
    return x, len(rows) - 1, rows


if __name__ == '__main__':
    dist = law()
    Ds = np.arange(2, DMAX + 1, 2)
    P = np.nonzero(C.primes_mask(10 ** 9, 10 ** 8))[0].astype(np.int64) + 10 ** 9
    edges = list(range(2, DMAX + 1, 12))                         # bands [lo, lo + 10]
    one = dist[1] / dist[1].sum()                                # class-averaged gap law (renormalised)
    res = {'lambda': LAM, 'bands': 'width 12', 'k': {}}
    for k in KS:
        d = P[k:] - P[:-k]; d = d[d <= DMAX]; N = len(d)
        cnt = np.array([np.sum((d >= lo) & (d <= lo + 10)) for lo in edges], float)
        lawk = dist[k] / dist[k].sum()
        e_law = np.array([N * lawk[(Ds >= lo) & (Ds <= lo + 10)].sum() for lo in edges])
        # Gamma(k) with mean k lam on even D (Gallagher), normalised within D <= DMAX
        gk = stats.gamma.pdf(Ds, k, scale=LAM); gk /= gk.sum()
        e_gam = np.array([N * gk[(Ds >= lo) & (Ds <= lo + 10)].sum() for lo in edges])
        # sum of k independent class-averaged gaps: conv index j corresponds to distance 2(j + k)
        if k == 1:
            ind = one.copy()
        else:
            full = one.copy()
            for _ in range(k - 1):
                full = np.convolve(full, one)
            dd = 2 * (np.arange(len(full)) + k); ind = np.array([full[dd == D].sum() for D in Ds]); ind /= ind.sum()
        e_ind = np.array([N * ind[(Ds >= lo) & (Ds <= lo + 10)].sum() for lo in edges])
        out = {'N': N}
        for name, e in (('law', e_law), ('gamma', e_gam), ('independent_gaps', e_ind)):
            x, dof, rows = chi2(cnt, e)
            out[name] = {'chi2': x, 'dof': dof, 'q999': float(stats.chi2.ppf(0.999, dof)),
                         'max_abs_z': float(max(abs(a - b) / math.sqrt(b) for a, b in rows))}
        out['bands'] = [{'lo': lo, 'count': int(c), 'law': float(e)} for lo, c, e in zip(edges, cnt, e_law)]
        res['k'][str(k)] = out
        print(k, 'N', N, {n: (round(out[n]['chi2'], 1), out[n]['dof'], round(out[n]['max_abs_z'], 1)) for n in ('law', 'gamma', 'independent_gaps')}, 'q999', round(out['law']['q999'], 1), flush=True)
    json.dump(res, open(os.path.join(OUT, 'r83.json'), 'w'), indent=1)
