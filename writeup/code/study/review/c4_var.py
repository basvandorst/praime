"""Check 4: Var/mean of prime counts in [L + jH, L + (j+1)H) under HL pairs, derived independently.

Count c = sum_{a<H} X_a. Conditional on n = L + jH = 0 (mod 10) (L and H multiples of 10):
  E X_a = S_0({a})/lam, E X_a X_b = S_0({a,b})/lam^2, with S_0 the singular series conditioned on n = 0 mod 10:
  local factors at p = 2, 5 replaced by 1[n+a, n+b units mod p]/(1-1/p)^k.  S_0 built from my own S (myS.py).
Var = sum_a EX_a + 2 sum_{a<b} E X_a X_b - (sum_a E X_a)^2 ; ratio Var/mean averaged over the window as
[<H/lam> + <(2s - H^2)/lam^2>] / <H/lam> with t uniform (every interval counts once) -- I check that this equals the
weight-1/ln t average of the ratio used in the write-up.
Also: the unaligned exact HL pair formula 1 - 1/lam + 2/(H lam) sum (H-h)(S(0,h)-1), and the M-S asymptotic.
"""
import math
import numpy as np
from myS import S

WIN = {'train': (10**7, 9 * 10**7), 'A': (10**10, 5 * 10**7), 'B': (10**12, 5 * 10**7), 'C': (10**15, 2 * 10**7)}


def S0_pair(a, b):
    # conditioned on n = 0 mod 10
    for p in (2, 5):
        if a % p == 0 or b % p == 0:
            return 0.0
    if a == b:
        raise ValueError
    full = S([a, b])
    # remove unconditional local factors at 2 and 5, insert conditional ones
    fac = 1.0
    for p in (2, 5):
        nu = len({a % p, b % p})
        fac *= (1.0 / (1 - 1 / p) ** 2) / ((1 - nu / p) / (1 - 1 / p) ** 2)
    return full * fac


def window_ratio(L, W, H, s, n=4001):
    t = np.linspace(L, L + W, n)
    lam = np.log(t)
    mean = H / lam
    var = H / lam + (2 * s - H * H) / lam ** 2
    r1 = np.trapezoid(var, t) / np.trapezoid(mean, t)              # pooled over intervals (uniform in t)
    r2 = np.trapezoid((1 / lam) * var / mean, t) / np.trapezoid(1 / lam, t)   # ratio averaged with weight 1/ln t
    return r1, r2


Spair = {}
for H in (50, 500):
    s = 0.0
    single = 0.0
    for a in range(H):
        if a % 2 and a % 5:
            single += 1 / (0.5 * 0.8)
    for h in range(1, H):
        if h % 2:
            continue
        if h not in Spair:
            Spair[h] = None
    # sum over pairs using the structure S_0(a, a+h) depends on h only (if both units): compute per h
    cache = {}
    for a in range(H):
        for b in range(a + 1, H):
            if a % 2 == 0 or b % 2 == 0 or a % 5 == 0 or b % 5 == 0:
                continue
            h = b - a
            key = (h, (a % 5, b % 5) if False else 0)
            if h not in cache:
                cache[h] = S0_pair(1, 1 + h) if (1 + h) % 5 else S0_pair(3, 3 + h)
            s += cache[h]
    # unaligned exact HL formula
    ps = sum((H - h) * (S([0, h]) - 1) if h % 2 == 0 else (H - h) * (0 - 1) for h in range(1, H))
    print('H=%d  sum_a E X_a * lam = %.6f (should be H)   s = %.6f' % (H, single, s))
    for w, (L, W) in WIN.items():
        r1, r2 = window_ratio(L, W, H, s)
        t = np.linspace(L, L + W, 4001); lam = np.log(t)
        un = np.trapezoid(1 / lam * (1 - 1 / lam + 2 * ps / (H * lam)), t) / np.trapezoid(1 / lam, t)
        ms = np.trapezoid(1 / lam * (1 - (math.log(H) - 1 + 0.5772156649 + math.log(2 * math.pi)) / lam), t) / np.trapezoid(1 / lam, t)
        print('  %-5s aligned (pooled) %.5f  aligned (1/ln-weighted ratio) %.5f   unaligned HL %.5f   M-S %.5f' % (w, r1, r2, un, ms))
