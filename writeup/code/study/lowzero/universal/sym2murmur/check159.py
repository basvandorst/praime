# Independent check of the data of round 159 at one holdout level: s(p) = sum_f lambda_f(p)^2 over the newforms of level 211,
# weight 2 and 4, against trace(T_p^2)/p^(k-1) from PARI's Hecke matrices on the new space (no trace formula, no class numbers).
import sys, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800_000_000)
nchk = nfail = 0; worst = 0.0
for k, N, fn in ((2, 211, 'data_k2_N211.npz'), (4, 173, 'data_k4_N173.npz')):
    z = np.load(sys.argv[1] + '/' + fn); P = z['p']; s = z['s']
    mf = pari.mfinit([N, k], 0)                      # 0: the new space
    assert int(pari.mfdim(mf)) == int(z['dim'])
    for p in (2, 3, 5, 7, 11, 13, 101, 211 if k == 4 else 199, 307, 401):
        if p == N: continue
        T = pari.mfheckemat(mf, p); tr = pari.trace(T * T)
        mine = float(tr) / p ** (k - 1); stored = float(s[list(P).index(p)])
        d = abs(mine - stored) / max(1.0, abs(mine)); worst = max(worst, d); nchk += 1; nfail += (d > 1e-9)
    print('k = %d, N = %d, dim %d: 10 primes compared' % (k, N, int(z['dim'])))
print('largest relative difference %.1e' % worst)
print('%d checks, %d failed' % (nchk, nfail))
