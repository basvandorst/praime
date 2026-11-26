"""Referee step 4: window means of t_1/dim and t_2/dim (N = 2q, k = 2) from PARI trace forms; noise level from the
per-prime scatter (SE of a log-weighted window mean); scaling with N; correlation between neighbouring levels."""
import math, time
import numpy as np
import cypari2
from sympy import primerange
pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6)
WIN = [(0.25 * i, 0.25 * (i + 1)) for i in range(16)]

def wstats(ps, N, X):
    y = ps / N; w = np.log(ps); m_ = []; se = []
    for a, b in WIN:
        m = (y > a) & (y <= b); ww = w[m]; x = X[m]; mu = np.sum(ww * x) / np.sum(ww)
        m_.append(mu); se.append(math.sqrt(np.sum(ww ** 2 * (x - mu) ** 2)) / np.sum(ww))
    return np.array(m_), np.array(se)

def main():
    groups = {'1e3': (499, 503, 509), '1e4': (5021, 5023, 5039), '2e4': (10103, 10111, 9973)}
    t0 = time.time()
    summ = {}
    for g, qs in groups.items():
        rows = []
        for q in qs:
            N = 2 * q
            ps = np.array([p for p in primerange(3, 4 * N + 1) if p != q], dtype=np.int64)
            dim = int(pari('mfdim([%d,2],0)' % N))
            tf = pari('mfcoefs(mftraceform([%d,2],0),%d)' % (N, 8 * N + 8))
            t1 = np.array([float(tf[int(p)]) for p in ps]) / dim; t2 = np.array([-float(tf[2 * int(p)]) for p in ps]) / dim
            r = {}
            for name, X in (('t1', t1), ('t2', t2), ('sym+', t1 + t2), ('sym-', t1 - t2)):
                r[name] = wstats(ps, N, X)
            rows.append((N, r))
        out = {}
        for name in ('t1', 't2', 'sym+', 'sym-'):
            M = np.array([r[name][0] for N, r in rows]); S = np.array([r[name][1] for N, r in rows])
            rms = float(np.sqrt(np.mean(M ** 2))); rse = float(np.sqrt(np.mean(S ** 2)))
            corr = float(np.mean([np.corrcoef(M[a], M[b])[0, 1] for a in range(3) for b in range(a + 1, 3)]))
            chi = float(np.sum((M / S) ** 2))
            out[name] = (rms, rse, corr, chi)
            # signal estimate: mean over the three levels vs its own noise
            mm = M.mean(axis=0); chi_mean = float(np.sum((mm / (np.sqrt(np.sum(S ** 2, axis=0)) / 3)) ** 2))
            out[name] += (chi_mean,)
        summ[g] = out
        print('%s levels %s (%.0fs)' % (g, [N for N, r in rows], time.time() - t0), flush=True)
        for name, v in out.items():
            print('   %-5s rms window mean %.4f | rms SE (noise) %.4f | mean corr between levels %+.2f | chi2(M/SE, 48 terms) %.1f | chi2 of 3-level mean (16) %.1f' % ((name,) + v))
    for name in ('t1', 'sym+', 'sym-'):
        a, b, c = (summ[g][name][0] for g in ('1e3', '1e4', '2e4'))
        print('%s: rms ratios 1e4/1e3 = %.3f (N^-1/2: %.3f), 2e4/1e4 = %.3f (N^-1/2: %.3f)' % (name, b / a, math.sqrt(1015 / 10055), c / b, math.sqrt(10055 / 20125)))

if __name__ == "__main__":
    main()
