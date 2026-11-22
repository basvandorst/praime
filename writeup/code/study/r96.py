"""Round 96 (preregister/round96.txt): R(X, delta) at 200 log-spaced X in [1e7, 1e10] from one pass of psi, and its
average over ln X, against the closed form D(delta) (round 95) and Hardy-Littlewood. usage: r96.py"""
import os, json, math, time
import numpy as np
from r91 import psi_at
import r93

OUT = r93.OUT
G0 = 0.5772156649015329; B = 1 - G0 - math.log(2 * math.pi)
D = lambda d: math.log(1 / d) - G0 - math.log(2 * math.pi) + 2 * math.log1p(d) + ((2 + d) / d) * math.log1p(d / 2)
HL = lambda d: math.log(1 / d) + B + (1 + d) * math.log1p(d) / d - 1
DELTAS = (0.1, 0.01)

if __name__ == '__main__':
    t0 = time.time(); N = 3_000_000
    u0, u1 = math.log(1e5), math.log(1e10); du = (u1 - u0) / N
    u = u0 + (np.arange(N) + np.random.default_rng(96).random(N)) * du
    x = np.exp(u); w = x * du
    q = np.concatenate([x] + [x * (1 + d) for d in DELTAS]); order = np.argsort(q, kind='stable')
    v = np.empty_like(q); v[order] = psi_at(q[order])
    ps0 = v[:N]
    Xs = np.exp(np.linspace(math.log(1e7), math.log(1e10), 200))
    res = {'points': N, 'X': Xs.tolist(), 'rows': {}}
    for k, d in enumerate(DELTAS):
        diff = v[(k + 1) * N:(k + 2) * N] - ps0 - d * x
        cum = np.cumsum(diff * diff * w)
        idx = np.searchsorted(x, Xs, side='right') - 1
        R = cum[idx] / (0.5 * d * Xs * Xs)
        # zero side for all X at once: sum_jk c_j conj(c_k) X^(2+i(g_j-g_k))/(2+i(g_j-g_k)) = X^2 v^T A conj(v), v_j = X^(i g_j)
        c = ((1 + d) ** r93.RHO - 1) / r93.RHO; G = r93.GAM
        A = (c[:, None] * np.conj(c)[None, :]) / (2 + 1j * (G[:, None] - G[None, :]))
        tail = (2 + d) * (math.log(r93.T / (2 * math.pi)) + 1) / (math.pi * r93.T * d)
        Rz = np.empty(len(Xs))
        for m, X in enumerate(Xs):
            vv = np.exp(1j * G * math.log(X))
            Rz[m] = 2 * float(np.real(vv @ (A @ np.conj(vv)))) / d + tail          # the -1 of the lower limit is negligible
        del A
        row = {'R_primes': R.tolist(), 'R_zeros': Rz.tolist(), 'avg_primes': float(np.mean(R)), 'avg_zeros': float(np.mean(Rz)),
               'D_closed': D(d), 'HL_exact': HL(d), 'max_abs_primes_minus_zeros': float(np.max(np.abs(R - Rz))),
               'min_R_primes': float(R.min()), 'max_R_primes': float(R.max())}
        res['rows'][str(d)] = row
        print('d %.2f  avg primes %.4f  avg zeros %.4f  D %.4f  HL %.4f  range %.3f..%.3f  max|p - z| %.4f' % (d, row['avg_primes'], row['avg_zeros'], D(d), HL(d), row['min_R_primes'], row['max_R_primes'], row['max_abs_primes_minus_zeros']), flush=True)
    res['seconds'] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, 'r96.json'), 'w'), indent=1)
