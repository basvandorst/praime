"""Round 90 (preregister/round90.txt): round 81's k-th neighbour law at small x, averaged over the window [1e7, 1e8)
(three lam, weights e^u/u), against the primes there; k = 1..4, q = 3..13, D <= 200. usage: r90.py"""
import os, json, math, time
import numpy as np
import common as C
import r47

assert r47.Y == 13
OUT = r47.OUT
NCLS = sum(1 for a in range(r47.PP) if math.gcd(a, r47.PP) == 1)
L0, L1 = 10 ** 7, 10 ** 8
LAMS = [math.log(L0), math.log(math.sqrt(L0 * L1)), math.log(L1)]
DMAX = 200; NZ = 48; KS = (1, 2, 3, 4); QS = (3, 5, 7, 11, 13)
Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)

if __name__ == '__main__':
    t0 = time.time(); dist = {k: np.zeros((len(LAMS), DMAX // 2)) for k in KS}
    for D in range(2, DMAX + 1, 2):
        st = r47.setup(D); M = st['M']; n = M.shape[1]
        for li, lam in enumerate(LAMS):
            rho0 = 1 / (lam * r47.RHO_DEN); p = rho0 * st['s1']
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
                dist[k][li, D // 2 - 1] = rho0 * st['S2'] * coef[:, k - 1].sum() / NCLS
    u = np.linspace(LAMS[0], LAMS[-1], 801); wt = np.exp(u) / u
    Ds = np.arange(2, DMAX + 1, 2)
    P = np.nonzero(C.primes_mask(L0, L1 - L0))[0].astype(np.int64) + L0
    res = {'lams': LAMS, 'DMAX': DMAX, 'primes': len(P), 'theory': {}, 'theory_mid_lam': {}, 'data': {}}
    for k in KS:
        avg = np.array([np.sum(np.interp(u, LAMS, dist[k][:, j]) * wt) / np.sum(wt) for j in range(DMAX // 2)])
        mid = dist[k][1]
        res['theory'][str(k)] = {str(q): float(avg[Ds % q == 0].sum() / avg.sum()) for q in QS}
        res['theory_mid_lam'][str(k)] = {str(q): float(mid[Ds % q == 0].sum() / mid.sum()) for q in QS}
        d = P[k:] - P[:-k]; d = d[d <= DMAX]; res['data'][str(k)] = {}
        for q in QS:
            m = float(np.mean(d % q == 0)); se = math.sqrt(m * (1 - m) / len(d))
            res['data'][str(k)][str(q)] = {'p': m, 'se': se, 'z': (res['theory'][str(k)][str(q)] - m) / se}
    zs = [v['z'] for k in res['data'] for v in res['data'][k].values()]
    res['max_abs_z'] = float(max(abs(z) for z in zs)); res['sum_z2'] = float(sum(z * z for z in zs))
    for k in KS:
        print(k, ' '.join('q%s %.4f/%.4f z%+.1f' % (q, res['theory'][str(k)][q], v['p'], v['z']) for q, v in res['data'][str(k)].items()))
    print('max |z| %.2f, sum z^2 %.1f for 20, %d s' % (res['max_abs_z'], res['sum_z2'], time.time() - t0), flush=True)
    json.dump(res, open(os.path.join(OUT, 'r90.json'), 'w'), indent=1)
