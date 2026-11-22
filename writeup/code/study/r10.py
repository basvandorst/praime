"""Round 10 (preregister/round10.txt): variance of prime counts in long intervals (new window V = [1e9, 2e9)) and
the third cumulant of counts in short intervals, both from exact Hardy-Littlewood sums on the aligned grid.
usage: r10.py var | skew"""
import sys, os, math, json, time
from itertools import combinations
import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import t as student_t
import common as C
import theory2 as T2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
C2 = 0.66016181584686957393


def aligned_s(H):
    """2 sum_{a<b} c R(b-a) u(a) u(b) / 2 = sum_h R(h) N_h c, by FFT (N_h) and a multiplicative loop over p (R)."""
    a = np.arange(H)
    u = (np.gcd(a, 10) == 1).astype(np.float64)
    Nh = np.rint(fftconvolve(u, u[::-1])[H:])            # lags 1..H-1
    R = np.full(H, C2 / (0.6 / 0.64))                    # index h (0 unused)
    for p in C.small_primes(H)[1:]:
        p = int(p)
        if p == 5:
            continue
        R[p::p] *= (p - 1) / (p - 2)
    return float(np.sum(Nh * R[1:]) / (0.5 * 0.8) ** 2)


def var_theory(H, L, W):
    s = aligned_s(H)
    return T2.window_avg(L, W, lambda lam: 1 + (2 * s - H * H) / (H * lam))


def var():
    L0, W0, K = 10 ** 9, 10 ** 9, 10
    Hs = [10 ** 4, 10 ** 5, 10 ** 6]
    assert abs(aligned_s(5000) - T2.aligned_sum(5000, 0)) < 1e-6 * abs(T2.aligned_sum(5000, 0))
    res = {str(H): {'meas': [], 'theory': []} for H in Hs}
    for b in range(K):
        Lb = L0 + b * (W0 // K)
        cnt = {H: [] for H in Hs}
        for j in range(10):                                # segments of 1e7 (a multiple of every H)
            M = C.primes_mask(Lb + j * 10 ** 7, 10 ** 7)
            for H in Hs[:2]:
                cnt[H].append(M.reshape(-1, H).sum(1))
            cnt[Hs[2]].append(M.reshape(-1, Hs[2]).sum(1))
        for H in Hs:
            c = np.concatenate(cnt[H]).astype(float)
            mid = Lb + np.arange(len(c)) * H + H / 2
            e = H / np.log(mid)
            res[str(H)]['meas'].append(float(np.var(c - e, ddof=1) / e.mean()))
            res[str(H)]['theory'].append(var_theory(H, Lb, W0 // K))
        print('block', b, {H: round(res[str(H)]['meas'][-1], 4) for H in Hs}, flush=True)
    for H in Hs:
        d = res[str(H)]
        r = np.array(d['meas']) - np.array(d['theory'])
        d['mean_meas'] = float(np.mean(d['meas'])); d['mean_theory'] = float(np.mean(d['theory']))
        d['t9'] = float(r.mean() / (r.std(ddof=1) / math.sqrt(K))); d['p'] = float(2 * student_t.sf(abs(d['t9']), K - 1))
        print(H, {k: (round(v, 4) if isinstance(v, float) else '') for k, v in d.items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'r10_var.json'), 'w'), indent=1)


def triple_sums(H):
    S0 = T2.S_cond(0, 10)
    U = np.array([a for a in range(H) if math.gcd(a, 10) == 1], dtype=np.int64)
    i, j = np.triu_indices(len(U), 1)
    s2 = float(S0(np.stack([U[i], U[j]], 1)).sum())
    s3 = 0.0
    it = combinations(range(len(U)), 3)
    while True:
        ch = np.array([x for _, x in zip(range(300000), it)], dtype=np.int64)
        if len(ch) == 0:
            break
        s3 += float(S0(U[ch]).sum())
    s1 = float(S0(U[:, None]).sum())
    return s1, s2, s3


def skew_theory(sums, L, W):
    s1, s2, s3 = sums

    def f(lam):
        m1 = s1 / lam
        m2 = m1 + 2 * s2 / lam ** 2
        m3 = m1 + 6 * s2 / lam ** 2 + 6 * s3 / lam ** 3
        return (m3 - 3 * m1 * m2 + 2 * m1 ** 3) / m1
    return T2.window_avg(L, W, f)


def skew():
    wins = dict(C.WINDOWS); wins['C'] = (10 ** 15, 2 * 10 ** 7)
    out = {}
    sums = {}
    for H in (50, 500):
        t0 = time.time(); sums[H] = triple_sums(H)
        print('sums', H, sums[H], round(time.time() - t0), flush=True)
    for win, (L, W) in wins.items():
        M = C.primes_mask(L, W)
        out[win] = {}
        for H in (50, 500):
            meas, th = [], []
            for b in range(10):
                Lb, Wb = L + b * (W // 10), W // 10
                c = M[Lb - L:Lb - L + Wb].reshape(-1, H).sum(1).astype(float)
                mid = Lb + np.arange(len(c)) * H + H / 2
                r = c - H / np.log(mid)
                meas.append(float(np.mean((r - r.mean()) ** 3) / np.mean(H / np.log(mid))))
                th.append(skew_theory(sums[H], Lb, Wb))
            rr = np.array(meas) - np.array(th)
            t9 = float(rr.mean() / (rr.std(ddof=1) / math.sqrt(10)))
            out[win][str(H)] = {'meas': float(np.mean(meas)), 'theory': float(np.mean(th)), 't9': t9,
                                'p': float(2 * student_t.sf(abs(t9), 9)),
                                'poisson': 1.0}
            print(win, H, out[win][str(H)], flush=True)
    out['_sums'] = {str(H): sums[H] for H in sums}
    json.dump(out, open(os.path.join(OUT, 'r10_skew.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'var': var, 'skew': skew}[sys.argv[1]]()
