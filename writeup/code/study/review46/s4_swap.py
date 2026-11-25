"""(B) swapped-term predictor variants for the near pairs (those with a frozen swap prediction).
mode 'asym' : all pairs, endpoint expansion (orders 2 and 3) with H, H' at 0, +-w from Euler products to 2e5 (and 5e4),
             plus H == 1 (product of the two one-point Gram terms) and H == H(0).
mode 'full' : a subset, full u-integral with the Euler product to 5e4 (du 0.06 = frozen setting, and du 0.02) and 2e5."""
import sys, time, math
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
from swapmine import *

mode = sys.argv[1]
S2 = np.load(OUT + '/s2_pairs.npz'); S1 = np.load(OUT + '/s1.npz'); Ds = list(S1['Ds'])
Xk = np.load(OUT + '/kron.npy'); P = np.load(OUT + '/primes.npy')
SW = S2['SW']; near = np.nonzero(np.any(SW != 0, axis=0))[0]
D1 = S2['D1'][near]; D2 = S2['D2'][near]; idx = {D: i for i, D in enumerate(Ds)}
print(mode, 'near pairs', len(near), flush=True)
t0 = time.time()
if mode == 'asym':
    out = {k: np.zeros((len(near), 16)) for k in ('a2_200k', 'a3_200k', 'a3_50k', 'h1', 'hconst')}
    H0s = np.zeros(len(near), complex); H0s50 = np.zeros(len(near), complex)
    for n, (a, b) in enumerate(zip(D1, D2)):
        a, b = int(a), int(b)
        psi = Xk[idx[a]].astype(int) * Xk[idx[b]].astype(int)
        kap = phi_ratio(abs(a)) * phi_ratio(abs(b))
        dL, setup = pair_setup(a, b, W)
        us = np.concatenate([[0.0], W, -W])
        for PM, tag in ((200000, '200k'), (50000, '50k')):
            m = P <= PM; Pf, s = hp_factors(psi[m], P[m])
            Hv, dH = H_and_dH(us, Pf, s)
            H0, dH0 = Hv[0], dH[0]; Hw, dHw = Hv[1:17], dH[1:17]; Hmw, dHmw = Hv[17:], dH[17:]
            if tag == '200k':
                out['a2_200k'][n] = cswap_asym(kap, dL, setup, W, H0, dH0, Hw, dHw, Hmw, dHmw, order=2)
                out['a3_200k'][n] = cswap_asym(kap, dL, setup, W, H0, dH0, Hw, dHw, Hmw, dHmw, order=3)
                one = np.ones(16, complex); zero = np.zeros(16, complex)
                out['h1'][n] = cswap_asym(kap, dL, setup, W, 1.0 + 0j, 0j, one, zero, one, zero, order=3)
                out['hconst'][n] = cswap_asym(kap, dL, setup, W, H0, 0j, H0 * one, zero, H0 * one, zero, order=3)
                H0s[n] = H0
            else:
                out['a3_50k'][n] = cswap_asym(kap, dL, setup, W, H0, dH0, Hw, dHw, Hmw, dHmw, order=3); H0s50[n] = H0
        if n % 500 == 0: print(n, '%.0fs' % (time.time() - t0), flush=True)
    np.savez(OUT + '/s4_asym.npz', near=near, D1=D1, D2=D2, H0=H0s, H0_50k=H0s50, **out)
    print('H(0) Euler 2e5 vs 5e4: rms rel diff %.4f' % np.sqrt(np.mean(np.abs(H0s - H0s50) ** 2) / np.mean(np.abs(H0s) ** 2)))
    print('done %.0fs' % (time.time() - t0))
else:
    nsub = int(sys.argv[2]); rng = np.random.default_rng(7)
    # subset: the pairs carrying the largest frozen predictions (half) and a random half
    big = np.argsort(-np.abs(SW[:, near]).max(0))[: nsub // 2]
    rest = np.setdiff1d(np.arange(len(near)), big); sub = np.sort(np.concatenate([big, rng.choice(rest, nsub - len(big), replace=False)]))
    out = {k: np.zeros((len(sub), 16)) for k in ('f50k_du02', 'f200k_du02', 'f10k_du02')}
    for i, n in enumerate(sub):
        a, b = int(D1[n]), int(D2[n])
        psi = Xk[idx[a]].astype(int) * Xk[idx[b]].astype(int)
        kap = phi_ratio(abs(a)) * phi_ratio(abs(b))
        dL, setup = pair_setup(a, b, W)
        for PM, du, tag in ((50000, 0.03, 'f50k_du02'), (200000, 0.03, 'f200k_du02'), (10000, 0.03, 'f10k_du02')):
            m = P <= PM; Pf, s = hp_factors(psi[m], P[m])
            out[tag][i] = cswap_full(kap, dL, setup, W, Pf, s, du=du)
        if i % 50 == 0: print(i, '%.0fs' % (time.time() - t0), flush=True)
    np.savez(OUT + '/s4_full_%d.npz' % nsub, sub=sub, near=near, **out)
    fr = SW[:, near[sub]].T
    for k, v in out.items():
        d = v - fr; print('%-11s vs frozen: rms diff %.5f, slope of frozen on it %.3f, slope of it on frozen %.3f' % (k, np.sqrt((d ** 2).mean()), (fr * v).sum() / (v * v).sum(), (fr * v).sum() / (fr * fr).sum()))
    print('done %.0fs' % (time.time() - t0))
