"""Round 81 (preregister/round81.txt): the k-th neighbour law. For class a and distance D, the number of primes strictly
between p and p + D has (order 2 of the Bahadur expansion around the conditioned base {0, D}, y = 13) the generating
function G_a(z) = exp( sum_i M_ai ln(1 + p_i (z-1)) + (1/2) sum_{i!=j} M_ai M_aj R_ij w_i w_j ), w_i = p_i (z-1)/(1 + p_i (z-1));
P_a(p + D is the k-th prime after p) = rho0 S2 [z^{k-1}] G_a(z), by a discrete Fourier sum over NZ points.
usage: r81.py"""
import os, json, math, time
import numpy as np
import common as C
import r47

assert r47.Y == 13
OUT = r47.OUT
PP = r47.PP
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], dtype=np.int64)
LAM = math.log(1.05e9); RHO0 = 1 / (LAM * r47.RHO_DEN)
DMAX = 240; NZ = 48; KS = (1, 2, 3, 4)
QS = (3, 5, 7, 11, 13)
Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)

if __name__ == '__main__':
    t0 = time.time()
    dist = {k: np.zeros(DMAX // 2) for k in KS}               # summed over classes, by D
    total_pair = np.zeros(DMAX // 2)
    for D in range(2, DMAX + 1, 2):
        st = r47.setup(D); M = st['M']; n = M.shape[1]; p = RHO0 * st['s1']
        logG = np.zeros((M.shape[0], NZ), complex)
        if n:
            lz = np.log1p(p[:, None] * (Z[None, :] - 1))          # (n, NZ)
            logG += M @ lz
        if n > 1:
            W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))   # (n, NZ)
            for m in range(NZ):
                X = M * W[:, m][None, :]                           # (classes, n) complex
                V = X.real @ st['R'] + 1j * (X.imag @ st['R'])
                logG[:, m] += 0.5 * np.sum(V * X, axis=1)
        G = np.exp(logG)
        coef = np.fft.fft(G, axis=1).real / NZ                    # c_j = (1/NZ) sum_m G(z_m) z_m^{-j}
        base = RHO0 * st['S2']
        for k in KS:
            dist[k][D // 2 - 1] = base * coef[:, k - 1].sum()
        total_pair[D // 2 - 1] = base * G[:, 0].real.sum()          # G(1) = 1: the pair density
        if D % 60 == 0:
            print('D', D, round(time.time() - t0), 's', flush=True)
    ncls = len(CLS)
    Ds = np.arange(2, DMAX + 1, 2)
    res = {'lambda': LAM, 'DMAX': DMAX, 'NZ': NZ, 'theory': {}, 'mass_within_DMAX': {}, 'pair_density_check': float(np.max(np.abs(total_pair / ncls - np.array([RHO0 * r47.S_gt(np.array([[0, d]]))[0] * sum(1 for a in CLS if math.gcd(int(a) + d, PP) == 1) / ncls for d in Ds]))))}
    for k in KS:
        dk = dist[k] / ncls
        res['mass_within_DMAX'][str(k)] = float(dk.sum())
        res['theory'][str(k)] = {str(q): float(dk[Ds % q == 0].sum() / dk.sum()) for q in QS}
    P = np.nonzero(C.primes_mask(10 ** 9, 10 ** 8))[0].astype(np.int64) + 10 ** 9
    res['data'] = {}
    for k in KS:
        d = P[k:] - P[:-k]; res.setdefault('data_mass_within_DMAX', {})[str(k)] = float(np.mean(d <= DMAX)); d = d[d <= DMAX]; res['data'][str(k)] = {}
        for q in QS:
            m = float(np.mean(d % q == 0)); se = math.sqrt(m * (1 - m) / len(d))
            res['data'][str(k)][str(q)] = {'p': m, 'se': se, 'z': (res['theory'][str(k)][str(q)] - m) / se}
    for k in KS:
        print(k, 'mass %.5f' % res['mass_within_DMAX'][str(k)], ' '.join('q%s %.4f/%.4f z%+.1f' % (q, res['theory'][str(k)][q], v['p'], v['z']) for q, v in res['data'][str(k)].items()), flush=True)
    print('pair density check', res['pair_density_check'])
    json.dump(res, open(os.path.join(OUT, 'r81.json'), 'w'), indent=1)
