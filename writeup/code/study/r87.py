"""Round 87 (preregister/round87.txt): round 81's k-th neighbour law (formula copied unchanged) at 1e18, against the
48,428 consecutive primes in [1e18, 1e18 + 2e6] (r17_primes.txt). usage: r87.py"""
import os, json, math, time
import numpy as np
import r47

assert r47.Y == 13
OUT = r47.OUT
NCLS = sum(1 for a in range(r47.PP) if math.gcd(a, r47.PP) == 1)
LAM = math.log(1e18 + 1e6); RHO0 = 1 / (LAM * r47.RHO_DEN)
DMAX = 480; NZ = 48; KS = tuple(range(1, 7)); QS = (3, 5, 7, 11, 13, 17, 19, 23)
Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)

if __name__ == '__main__':
    t0 = time.time(); dist = {k: np.zeros(DMAX // 2) for k in KS}
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
        if D % 120 == 0:
            print('D', D, round(time.time() - t0), 's', flush=True)
    Ds = np.arange(2, DMAX + 1, 2)
    P = np.array([int(x) for x in open(os.path.join(OUT, 'r17_primes.txt'))], dtype=np.int64)
    assert np.all(np.diff(P) > 0)
    res = {'lambda': LAM, 'DMAX': DMAX, 'primes': len(P), 'theory': {}, 'data': {}, 'mass_within_DMAX': {}}
    for k in KS:
        dk = dist[k]; res['mass_within_DMAX'][str(k)] = float(dk.sum())
        res['theory'][str(k)] = {str(q): float(dk[Ds % q == 0].sum() / dk.sum()) for q in QS}
        d = P[k:] - P[:-k]; d = d[d <= DMAX]; res['data'][str(k)] = {}
        for q in QS:
            m = float(np.mean(d % q == 0)); se = math.sqrt(m * (1 - m) / len(d))
            res['data'][str(k)][str(q)] = {'p': m, 'se': se, 'z': (res['theory'][str(k)][str(q)] - m) / se,
                                           'z_uniform': (1 / (q - 1) - m) / se}
    zs = [v['z'] for k in res['data'] for v in res['data'][k].values()]
    res['max_abs_z'] = float(max(abs(z) for z in zs)); res['sum_z2'] = float(sum(z * z for z in zs))
    for k in KS:
        print(k, ' '.join('q%s %.4f/%.4f z%+.1f u%+.1f' % (q, res['theory'][str(k)][q], v['p'], v['z'], v['z_uniform']) for q, v in res['data'][str(k)].items()))
    print('max |z| %.2f, sum z^2 %.1f for 48' % (res['max_abs_z'], res['sum_z2']), flush=True)
    json.dump(res, open(os.path.join(OUT, 'r87.json'), 'w'), indent=1)
