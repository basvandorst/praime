"""Round 42 (preregister/round42.txt): the spread of twin counts across the admissible classes mod Q (Q = 210: 15
classes, Q = 2310: 135 classes), all with the same law, against the Hardy-Littlewood four-point sums. Calibration-free:
the measured statistic is Var_a(X_a)/mean over the classes (ddof 1). usage: r42.py"""
import os, json, math
import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import chi2
import common as C
import theory_r3 as R3

PR = C.small_primes(20_000_000).astype(np.int64)
WINS = {'train': (10 ** 7, 9 * 10 ** 7), '1e8-1e9': (10 ** 8, 9 * 10 ** 8), 'A': (10 ** 10, 5 * 10 ** 7),
        'B': (10 ** 12, 5 * 10 ** 7), 'C': (10 ** 15, 2 * 10 ** 7)}


def rho_twin(n, Q):
    """rho(k) = S({0,2,Qk,Qk+2}) / S({0,2})^2 with the factors p | Q removed (they cancel), k = 0..n-1."""
    qf = [int(p) for p in PR[PR < 50] if Q % p == 0]
    P = PR[PR > max(qf)]
    pf = P.astype(np.float64)
    fn = lambda p: (1 - 4 / p) / (1 - 2 / p) ** 2
    L = np.full(n, float(np.sum(np.log(fn(pf)))))
    for p in P[P < n]:
        p = int(p)
        L[p::p] += math.log((1 / (1 - 2 / p)) / fn(p))
    k = np.arange(n, dtype=np.int64)
    lim = int(math.isqrt(Q * n + 2)) + 1
    for s in (-1, 1):
        cof = np.abs(Q * k + 2 * s)
        for _ in range(64):                                   # Qk +- 2 is even; its odd part has no factor of Q
            m = (cof % 2 == 0) & (cof > 0)
            if not m.any():
                break
            cof[m] //= 2
        for p in P[P <= lim]:
            p = int(p)
            k0 = (-2 * s * pow(Q, -1, p)) % p
            idx = np.arange(k0, n, p)
            L[idx] += math.log(((1 - 3 / p) / (1 - 2 / p) ** 2) / fn(p))
            for _ in range(40):
                m = cof[idx] % p == 0
                if not m.any():
                    break
                cof[idx[m]] //= p
        big = cof > 1
        pb = cof[big].astype(np.float64)
        L[big] += np.log(((1 - 3 / pb) / (1 - 2 / pb) ** 2) / fn(pb))
    return np.exp(L)


def predict(L0, W, Q):
    qf = [int(p) for p in PR[PR < 50] if Q % p == 0]
    P = PR[PR > max(qf)].astype(np.float64)
    c = math.prod(1 / (1 - 1 / p) ** 2 for p in qf) * math.exp(float(np.sum(np.log((1 - 2 / P) / (1 - 1 / P) ** 2))))
    n = W // Q
    m = L0 + Q * np.arange(n) + Q / 2
    r = c / (np.log(m) * np.log(m + 2))
    w = fftconvolve(r, r[::-1])[n - 1:]
    rho = rho_twin(n, Q)
    var = float(np.sum(r * (1 - r))) + 2 * float(np.sum(w[1:] * (rho[1:] - 1)))
    return var / float(r.sum())


def count(L0, W, Qs):
    out = {Q: np.zeros(Q) for Q in Qs}
    last = None; seg = 10 ** 7
    for s in range(0, W, seg):
        P = np.nonzero(C.primes_mask(L0 + s, min(seg, W - s)))[0].astype(np.int64) + L0 + s
        if last is not None:
            P = np.concatenate([[last], P])
        last = int(P[-1])
        t = P[:-1][np.diff(P) == 2]
        for Q in Qs:
            out[Q] += np.bincount(t % Q, minlength=Q)
    return out


if __name__ == '__main__':
    res = {}
    for win, (L0, W) in WINS.items():
        cnt = count(L0, W, (210, 2310))
        res[win] = {}
        for Q in (210, 2310):
            adm = [a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1]
            X = cnt[Q][adm]
            K = len(adm)
            meas = float(np.var(X, ddof=1) / X.mean())
            pred = predict(L0, W, Q)
            # lower-tail probability of the measured spread if K-1 times spread/pred is chi^2 with K-1 dof
            p_low = float(chi2.cdf((K - 1) * meas / pred, K - 1))
            res[win][str(Q)] = {'classes': K, 'twins': int(X.sum()), 'measured': meas, 'predicted': pred,
                                'ratio': meas / pred, 'se_ratio': math.sqrt(2 / (K - 1)), 'p_lower': p_low}
            print(win, Q, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[win][str(Q)].items()}, flush=True)
    for Q in ('210', '2310'):
        rs = np.array([res[w][Q]['ratio'] for w in WINS]); se = np.array([res[w][Q]['se_ratio'] for w in WINS])
        wt = 1 / se ** 2
        res['pooled_' + Q] = {'ratio': float(np.sum(wt * rs) / wt.sum()), 'se': float(1 / math.sqrt(wt.sum())),
                              'fisher_p': float(chi2.sf(-2 * sum(math.log(res[w][Q]['p_lower']) for w in WINS), 2 * len(WINS)))}
        print('pooled', Q, res['pooled_' + Q], flush=True)
    json.dump(res, open(os.path.join(R3.OUT, 'r42.json'), 'w'), indent=1)
