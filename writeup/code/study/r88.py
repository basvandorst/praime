"""Round 88 (preregister/round88.txt): round 81's k-th neighbour law recomputed with y = 17 (classes mod 510510),
compared with round 81 (y = 13) and with the data of r81.json. usage: r88.py"""
import os, json, math, time
os.environ['R47_Y'] = '17'
import numpy as np
import r47

assert r47.Y == 17
OUT = r47.OUT
NCLS = sum(1 for a in range(r47.PP) if math.gcd(a, r47.PP) == 1)
LAM = math.log(1.05e9); RHO0 = 1 / (LAM * r47.RHO_DEN)
DMAX = 240; NZ = 48; KS = (1, 2, 3, 4); QS = (3, 5, 7, 11, 13)
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
        if D % 60 == 0:
            print('D', D, round(time.time() - t0), 's', flush=True)
    Ds = np.arange(2, DMAX + 1, 2); r81 = json.load(open(os.path.join(OUT, 'r81.json')))
    res = {'y': 17, 'lambda': LAM, 'theory': {}, 'diff_over_se': {}, 'z_data': {}, 'mass_within_DMAX': {}}
    for k in KS:
        dk = dist[k]; res['mass_within_DMAX'][str(k)] = float(dk.sum())
        res['theory'][str(k)] = {str(q): float(dk[Ds % q == 0].sum() / dk.sum()) for q in QS}
        res['diff_over_se'][str(k)] = {str(q): (res['theory'][str(k)][str(q)] - r81['theory'][str(k)][str(q)]) / r81['data'][str(k)][str(q)]['se'] for q in QS}
        res['z_data'][str(k)] = {str(q): (res['theory'][str(k)][str(q)] - r81['data'][str(k)][str(q)]['p']) / r81['data'][str(k)][str(q)]['se'] for q in QS}
    res['max_abs_diff_over_se'] = max(abs(v) for d in res['diff_over_se'].values() for v in d.values())
    res['max_abs_z_data'] = max(abs(v) for d in res['z_data'].values() for v in d.values())
    for k in KS:
        print(k, ' '.join('q%s %.5f d/se %+.2f z %+.1f' % (q, res['theory'][str(k)][q], res['diff_over_se'][str(k)][q], res['z_data'][str(k)][q]) for q in map(str, QS)))
    print('max |diff|/se %.3f, max |z| %.2f, %d s' % (res['max_abs_diff_over_se'], res['max_abs_z_data'], time.time() - t0), flush=True)
    json.dump(res, open(os.path.join(OUT, 'r88.json'), 'w'), indent=1)
