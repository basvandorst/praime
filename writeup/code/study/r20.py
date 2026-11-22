"""Round 20 (preregister/round20.txt): the distribution of the number of primes in [n, n + H) from the Bahadur form
of the probability generating function G(z) = E prod (1 - (1 - z) X_m), evaluated at roots of unity, P(k) by FFT."""
import os, sys, json, math, time
from itertools import combinations, permutations
import numpy as np
from scipy.stats import t as student_t
import theory2 as T2
import theory_r3 as R3
import common as C

OUT = R3.OUT
NZ = 64
KMAX = {30: 7, 60: 10}


class Interval:
    def __init__(self, H):
        self.H = H
        pts = np.arange(H, dtype=np.int64)
        i, j = np.triu_indices(H, 1)
        self.pairs = (i, j, T2.S_many(np.stack([pts[i], pts[j]], 1)))
        T = np.array(list(combinations(range(H), 3)), dtype=np.int64)
        s3v = T2.S_many(T)
        self.s3 = np.zeros((H, H, H))
        for pm in permutations(range(3)):
            self.s3[T[:, pm[0]], T[:, pm[1]], T[:, pm[2]]] = s3v
        self.trip = (T[:, 0], T[:, 1], T[:, 2], s3v)
        self.s2 = np.zeros((H, H)); self.s2[i, j] = self.pairs[2]; self.s2[j, i] = self.pairs[2]
        Q = np.array(list(combinations(range(H), 4)), dtype=np.int64)
        self.quad = (Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3], T2.S_many(Q))

    def logG(self, z, lam):
        w = 1 - z
        p = w / lam
        H = self.H
        out = H * np.log(1 - p)
        q = 1 - p
        i, j, s2v = self.pairs
        out += np.sum(w * w * s2v / lam ** 2 - p * p) / q ** 2
        a, b, c, s3v = self.trip
        m2 = w * w * self.s2 / lam ** 2
        k3 = w ** 3 * s3v / lam ** 3 - (m2[a, b] + m2[a, c] + m2[b, c]) * p + 2 * p ** 3
        out -= np.sum(k3) / q ** 3
        a, b, c, d, s4v = self.quad
        m3 = lambda x, y, zz: w ** 3 * self.s3[x, y, zz] / lam ** 3
        k4 = (w ** 4 * s4v / lam ** 4
              - (m3(a, b, c) + m3(a, b, d) + m3(a, c, d) + m3(b, c, d)) * p
              - m2[a, b] * m2[c, d] - m2[a, c] * m2[b, d] - m2[a, d] * m2[b, c]
              + 2 * (m2[a, b] + m2[a, c] + m2[a, d] + m2[b, c] + m2[b, d] + m2[c, d]) * p * p
              - 6 * p ** 4)
        out += np.sum(k4) / q ** 4
        return out

    def pk(self, lam):
        zs = np.exp(2j * np.pi * np.arange(NZ) / NZ)
        G = np.array([np.exp(self.logG(z, lam)) for z in zs])
        return np.real(np.fft.fft(G)) / NZ          # P(k) = (1/N) sum_j G(z_j) z_j^{-k}


def law(iv, L, W):
    u = np.linspace(math.log(L), math.log(L + W), 7)
    Pu = np.array([iv.pk(x) for x in u])
    wts = np.exp(u)                                  # every n counts: weight dt = e^u du
    return np.trapezoid(wts[:, None] * Pu, u, axis=0) / np.trapezoid(wts, u)


def counts(M, H):
    c = np.concatenate((np.zeros(1, dtype=np.int32), np.cumsum(M, dtype=np.int32)))
    return c[H:] - c[:-H]


def main():
    out = {}
    for H in [int(x) for x in os.environ.get('R20_H', '30,60').split(',')]:
        t0 = time.time(); iv = Interval(H); print('H', H, 'built', round(time.time() - t0), flush=True)
        for win, (L, W) in [(k, v) for k, v in R3.WIN.items() if k in os.environ.get('R20_WIN', 'train,A,B,C').split(',')]:
            M = C.primes_mask(L, W + H)
            K = KMAX[H]
            emp, th = [], []
            for b in range(10):                       # equal blocks; the whole window is their mean
                lo, n = b * W // 10, W // 10
                cnt = counts(M[lo:lo + n + H - 1], H)
                e = np.bincount(cnt, minlength=K + 1)[:K + 1] / len(cnt)
                emp.append(e); th.append(law(iv, L + lo, n)[:K + 1])
            emp = np.array(emp); th = np.array(th)
            emp = np.vstack([emp.mean(0), emp]); th = np.vstack([th.mean(0), th])
            lam = math.log(L + W / 2)
            mu = H / lam
            pois = np.array([math.exp(-mu) * mu ** k / math.factorial(k) for k in range(K + 1)])
            r = emp[1:] - th[1:]
            t9 = r.mean(0) / (r.std(0, ddof=1) / math.sqrt(10))
            res = {'H': H, 'data': emp[0].tolist(), 'law': th[0].tolist(), 'poisson': pois.tolist(), 't9': t9.tolist(),
                   'p': (2 * student_t.sf(np.abs(t9), 9)).tolist(), 'law_sum': float(np.sum(law(iv, L, W)))}
            out['%s_H%d' % (win, H)] = res
            print(win, H, 'k:', ' '.join('%d:%.4f/%.4f/%.4f(t%.1f)' % (k, emp[0][k], th[0][k], pois[k], t9[k]) for k in range(K + 1)), flush=True)
    json.dump(out, open(os.path.join(OUT, 'r20_H%s%s.json' % (os.environ.get('R20_H', 'all').replace(',', '_'), ('_' + os.environ['R20_WIN'].replace(',', '_')) if os.environ.get('R20_WIN') else '')), 'w'), indent=1)


if __name__ == '__main__':
    main()
