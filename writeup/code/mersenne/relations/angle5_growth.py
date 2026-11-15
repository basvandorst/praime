"""Angle 5: log2 a_n against n, and the point processes x = log2 a, y = log2 log2 a.

Residuals e_n of the least-squares line x_n ~ n:
  res_sd, res_ac1..3, res_runs (sign runs), res_g (Fisher's g of the periodogram)
Gaps g_n = x_{n+1} - x_n (51 gaps; Poisson in x <=> iid exponential gaps):
  gap_cv (CV > 1 clustering, < 1 repulsion), gap_cv2 (Holt's local CV2), gap_ac1,
  gap_lb5 (Ljung-Box Q(5)), gap_g (Fisher's g), gap_small (# gaps < 0.1 mean),
  gap_max (max gap / mean)
Point process y_n = log2 log2 a_n, n = 2..52 (inhomogeneous; the nulls calibrate):
  y_cv, y_cv2, y_pairs01/03 (# pairs with |y_i - y_j| < 0.01 / 0.03)
All two-sided except Fisher's g and Ljung-Box (hi).
"""
import json, sys, time
import numpy as np


def acf(x, lag):
    x = x - x.mean()
    return float((x[:-lag] * x[lag:]).sum() / (x * x).sum())


def fisher_g(x):
    x = x - x.mean()
    P = np.abs(np.fft.rfft(x)[1:]) ** 2
    if len(x) % 2 == 0:
        P = P[:-1]
    return float(P.max() / P.sum())


def cv2(g):
    return float(np.mean(2 * np.abs(g[1:] - g[:-1]) / (g[1:] + g[:-1])))


def all_stats(a):
    a = np.array(a, dtype=np.float64)
    x = np.log2(a)
    n = np.arange(1, len(a) + 1, dtype=float)
    c = np.polyfit(n, x, 1)
    e = x - np.polyval(c, n)
    out = {'res_sd': float(e.std()), 'res_ac1': acf(e, 1), 'res_ac2': acf(e, 2), 'res_ac3': acf(e, 3),
           'res_runs': int(1 + np.sum(np.sign(e[1:]) != np.sign(e[:-1]))), 'res_g': fisher_g(e)}
    g = np.diff(x)
    m = g.mean()
    out['gap_cv'] = float(g.std() / m)
    out['gap_cv2'] = cv2(g)
    out['gap_ac1'] = acf(g, 1)
    N = len(g)
    out['gap_lb5'] = float(N * (N + 2) * sum(acf(g, k) ** 2 / (N - k) for k in range(1, 6)))
    out['gap_g'] = fisher_g(g)
    out['gap_small'] = int(np.sum(g < 0.1 * m))
    out['gap_max'] = float(g.max() / m)
    y = np.log2(x[1:])
    h = np.diff(y)
    out['y_cv'] = float(h.std() / h.mean())
    out['y_cv2'] = cv2(h)
    dif = np.abs(y[:, None] - y[None, :])[np.triu_indices(len(y), 1)]
    out['y_pairs01'] = int(np.sum(dif < 0.01))
    out['y_pairs03'] = int(np.sum(dif < 0.03))
    return out


SIDES = {k: 'two' for k in ['res_sd', 'res_ac1', 'res_ac2', 'res_ac3', 'res_runs', 'gap_cv', 'gap_cv2',
                            'gap_ac1', 'gap_small', 'gap_max', 'y_cv', 'y_cv2', 'y_pairs01', 'y_pairs03']}
SIDES.update({'res_g': 'hi', 'gap_lb5': 'hi', 'gap_g': 'hi'})

if __name__ == '__main__':
    D = json.load(open('nulls.json'))
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    t0 = time.time()
    A = D['real']
    real = all_stats(A)
    with open('angle5_real.txt', 'w') as rep:
        rep.write('real stats: %s\n' % real)
        x = np.log2(np.array(A, float))
        rep.write('slope of log2 a_n vs n: %.4f bits/term (rate %.4f terms/bit; e^gamma = %.4f)\n'
                  % (np.polyfit(np.arange(1, 53), x, 1)[0], 1 / np.polyfit(np.arange(1, 53), x, 1)[0],
                     np.exp(np.euler_gamma)))
        rep.write('gaps in bits: %s\n' % np.round(np.diff(x), 4).tolist())
    out = {s: {'side': SIDES[s], 'real': real[s], 'a': [], 'b': []} for s in SIDES}
    for fam in ('a', 'b'):
        for seq in D[fam][:nmax]:
            st = all_stats(seq)
            for s in SIDES:
                out[s][fam].append(st[s])
    json.dump(out, open('angle5.json', 'w'))
    print('done', time.time() - t0)
