"""d = -239: F_3(1/3) = 0 (3 splits completely, lambda_3 = 2). Recompute c' = c/F_3(1/3) analytically (p = 3 factor of F removed),
the Artin lines at r = 3, 1/3 (-eps c' S_3(+-1)) and the M1 mechanism (-sum_k c rho(1/k), only 3 | k survive), against the frozen values."""
import sys, os, json, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from bz_lib import build
from ulib import local_poly, local_data, series_inv, primes_upto, to_c
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
pol = 'x^3 - x - 3'; d = -239
Lr = build(pari, ['artin', pol]); Lc = pari.lfuncreate(d)
A1 = to_c(pari.lfun(Lc, 1)).real * to_c(pari.lfun(Lr, 1)).real
lb = 0.0; data = {}
for p in primes_upto(3000):
    p = int(p); Pp = local_poly(pari, Lr, p, 1); F, S, lam = local_data(Pp, p); data[p] = (F, S, Pp)
    a = np.convolve(local_poly(pari, Lc, p, 1), Pp); Aval = 1 / sum(a[j] * p ** (-j) for j in range(len(a)))
    lb += math.log(abs(Aval / (1 - 1 / p))) + (0.0 if p == 3 else math.log(abs(F)))
    sgn = 1
cprime = math.exp(lb) / A1 * np.sign(np.prod([data[p][0].real for p in data if p != 3]))
F3, S3, P3 = data[3]
print('P_3 (analytic local polynomial) =', [complex(round(x.real, 6), round(x.imag, 6)) for x in P3], ' F_3(1/3) =', F3)
fr = json.load(open('constants/T5_d239.json'))
for v, nm in ((1, '3'), (-1, '1/3')):
    amp = -cprime * S3[v]
    print('line r=%-3s analytic amplitude per unit %+.6f%+.6fi   frozen %s' % (nm, amp.real, amp.imag, fr['rho']['lines'][nm]))
# M1 mechanism: -sum_{k<=K} c rho(1/k) = -sum_{k<=K, 3|k} c' * [p^{v/2} S_3(-v)] * prod_{q != 3} rho_q(1/q^w)
K = fr['mech']['K']; rho = {}
for p in primes_upto(K):
    p = int(p); F, S, Pp = data[p]; lam = series_inv(Pp, 80)
    for vv in range(1, 8):
        if p ** vv > K: break
        s = sum(np.conj(lam[h - vv]) * Pp[h] * p ** (-(2 * h - vv) / 2) for h in range(vv, len(Pp)))
        rho[(p, vv)] = (p ** (vv / 2) * s / F) if p != 3 else p ** (vv / 2) * s      # for p = 3: without 1/F_3 (absorbed in c')
M1 = 0.0
for k in range(1, K + 1):
    if k % 3: continue
    f = pari.factor(k); val = 1
    for i in range(len(f[0])): val *= rho[(int(f[0][i]), int(f[1][i]))]
    M1 += val.real
print('M1 mechanism analytic %.4f   frozen %.4f' % (-cprime * M1, fr['mech']['M1']))
