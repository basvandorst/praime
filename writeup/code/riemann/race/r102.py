"""Round 102 (preregister/round102.txt): R_chi(X, d) for the prime race (chi_4, chi_3) at 200 log-spaced X in [1e6, 1e9],
averaged over ln X, against the Guinand-Weil value D_chi(d) (dchi.py). usage: r102.py"""
import os, sys, json, math, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'study'))
import common as C

XMAX = 10 ** 9; QMAX = int(1.1 * XMAX) + 2
CHARS = {'chi4': (4, [0, 1, 0, -1]), 'chi3': (3, [0, 1, -1])}
DS = (0.1, 0.01)


def psi_chi_at(qs, q, chi):
    """psi(x, chi) at sorted points qs (all <= QMAX)."""
    out = np.zeros(len(qs)); base = 0.0; SEG = 10 ** 7
    pw = []
    for p in C.small_primes(int(math.isqrt(QMAX)) + 2).tolist():
        v = p * p; k = 2
        while v <= QMAX:
            pw.append((v, math.log(p) * chi[p % q] ** k)); v *= p; k += 1
    pw.sort(); pwv = np.array([a for a, _ in pw], float); pwc = np.cumsum([b for _, b in pw])
    for s0 in range(0, QMAX, SEG):
        n = min(SEG, QMAX - s0)
        P = np.nonzero(C.primes_mask(max(s0, 2), n - (2 if s0 == 0 else 0)))[0].astype(np.int64) + max(s0, 2)
        w = np.log(P.astype(float)) * np.array(chi, float)[P % q]
        cs = np.concatenate([[0.0], np.cumsum(w)])
        lo = np.searchsorted(qs, s0); hi = np.searchsorted(qs, s0 + n)
        if hi > lo:
            k = np.searchsorted(P, np.floor(qs[lo:hi]).astype(np.int64), side='right')
            out[lo:hi] = base + cs[k]
        base += float(np.sum(w))
    k = np.searchsorted(pwv, np.floor(qs), side='right')
    out += np.concatenate([[0.0], pwc])[k]
    return out


if __name__ == '__main__':
    t0 = time.time(); N = 2_000_000
    u0, u1 = math.log(1e4), math.log(XMAX); du = (u1 - u0) / N
    u = u0 + (np.arange(N) + np.random.default_rng(102).random(N)) * du
    x = np.exp(u); w = x * du
    Xs = np.exp(np.linspace(math.log(1e6), math.log(XMAX), 200))
    dch = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dchi.json')))
    res = {'points': N, 'rows': {}}
    for name, (q, chi) in CHARS.items():
        qq = np.concatenate([x] + [x * (1 + d) for d in DS]); order = np.argsort(qq, kind='stable')
        v = np.empty_like(qq); v[order] = psi_chi_at(qq[order], q, chi)
        p0 = v[:N]
        for kk, d in enumerate(DS):
            diff = v[(kk + 1) * N:(kk + 2) * N] - p0
            cum = np.cumsum(diff * diff * w); idx = np.searchsorted(x, Xs, side='right') - 1
            R = cum[idx] / (0.5 * d * Xs * Xs)
            Dv = dch[str(d)][name]
            res['rows']['%s_%s' % (name, d)] = {'avg': float(np.mean(R)), 'D': Dv, 'diff': float(np.mean(R)) - Dv,
                                                'with_pole': Dv + 2 * math.log1p(d), 'min': float(R.min()), 'max': float(R.max()), 'R': R.tolist()}
            print(name, d, 'avg %.4f  D %.4f  diff %+.4f  (with pole %.4f)  range %.3f..%.3f' % (np.mean(R), Dv, np.mean(R) - Dv, Dv + 2 * math.log1p(d), R.min(), R.max()), flush=True)
    res['seconds'] = time.time() - t0
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r102.json'), 'w'), indent=1)
