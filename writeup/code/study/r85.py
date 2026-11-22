"""Round 85 (preregister/round85.txt): the last-digit table of consecutive primes from round 81's k-th neighbour law,
split by the class of p mod 10: P(last digit of p_{n+k} = b | last digit of p_n = a), a, b in {1, 3, 7, 9}, k = 1..3,
conditioned on p_{n+k} - p_n <= DMAX, against the primes. usage: R85_L0=.. R85_L1=.. R85_DMAX=.. R85_OUT=.. r85.py"""
import os, json, math, time
import numpy as np
import common as C
import r47

assert r47.Y == 13
OUT = r47.OUT
PP = r47.PP
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], dtype=np.int64)
L0 = int(float(os.environ.get('R85_L0', 1e9))); L1 = int(float(os.environ.get('R85_L1', 1.1e9)))
LAM = math.log((L0 + L1) / 2); RHO0 = 1 / (LAM * r47.RHO_DEN)
DMAX = int(os.environ.get('R85_DMAX', 240)); NZ = 48; KS = (1, 2, 3)
DIG = (1, 3, 7, 9)
Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)
PATH = os.path.join(OUT, os.environ.get('R85_OUT', 'r85_1e9.json'))

if __name__ == '__main__':
    t0 = time.time()
    T = {k: np.zeros((4, 4)) for k in KS}                      # T[k][a, b]: summed law weight, a = digit of p, b = digit of p + D
    for D in range(2, DMAX + 1, 2):
        st = r47.setup(D); M = st['M']; n = M.shape[1]; p = RHO0 * st['s1']
        cls = CLS[np.gcd(CLS + D, PP) == 1]                       # the rows of M, in setup's order
        logG = np.zeros((M.shape[0], NZ), complex)
        if n:
            logG += M @ np.log1p(p[:, None] * (Z[None, :] - 1))
        if n > 1:
            W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))
            for m in range(NZ):
                X = M * W[:, m][None, :]
                logG[:, m] += 0.5 * np.sum((X.real @ st['R'] + 1j * (X.imag @ st['R'])) * X, axis=1)
        coef = np.fft.fft(np.exp(logG), axis=1).real / NZ
        da = np.array([DIG.index(x) for x in (cls % 10)]); db = np.array([DIG.index(x) for x in ((cls + D) % 10)])
        for k in KS:
            np.add.at(T[k], (da, db), RHO0 * st['S2'] * coef[:, k - 1])
    cnt = {k: np.zeros((4, 4)) for k in KS}; carry = np.zeros(0, np.int64)
    for s0 in range(L0, L1, 10 ** 8):
        P = np.concatenate([carry, np.nonzero(C.primes_mask(s0, min(10 ** 8, L1 - s0)))[0].astype(np.int64) + s0])
        for k in KS:
            a = P[:-k]; b = P[k:]
            if len(carry):
                a = a[max(0, len(carry) - k):]; b = b[max(0, len(carry) - k):]
            ok = (b - a) <= DMAX
            ia = np.searchsorted(DIG, a[ok] % 10); ib = np.searchsorted(DIG, b[ok] % 10)
            np.add.at(cnt[k], (ia, ib), 1)
        carry = P[-max(KS):]
    res = {'lambda': LAM, 'L0': L0, 'L1': L1, 'DMAX': DMAX, 'digits': DIG, 'k': {}}
    zall = []
    for k in KS:
        law = T[k] / T[k].sum(1, keepdims=True); N = cnt[k].sum(1, keepdims=True); dat = cnt[k] / N
        se = np.sqrt(dat * (1 - dat) / N); z = (law - dat) / se; zall += list(z.ravel())
        res['k'][str(k)] = {'law': law.tolist(), 'data': dat.tolist(), 'z': z.tolist(), 'row_counts': N.ravel().tolist()}
        print('k', k, 'max|z| %.2f' % np.max(np.abs(z)))
        for i, a in enumerate(DIG):
            print('  ', a, ' '.join('%.4f/%.4f(%+.1f)' % (law[i, j], dat[i, j], z[i, j]) for j in range(4)))
    res['max_abs_z'] = float(np.max(np.abs(zall))); res['chi2'] = float(np.sum(np.square(zall)))
    print('max |z| %.2f, chi2 %.1f over %d entries, %d s' % (res['max_abs_z'], res['chi2'], len(zall), time.time() - t0), flush=True)
    json.dump(res, open(PATH, 'w'), indent=1)
