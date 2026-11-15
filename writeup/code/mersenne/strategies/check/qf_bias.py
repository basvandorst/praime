"""Chebyshev bias in Chebotarev form for x^2 + P y^2 (P = 1..100): the primes represented (split completely in the ring
class field, degree 2h) fall short of pi/(2h). Prediction (Rubinstein-Sarnak style, prime squares): count ~ [Li(N) - (1/2)*s*Li(sqrt N)]/(2h)
with s = #{g in G : g^2 = 1} = h + a (a = ambiguous classes = |Cl[2]|), while pi(N) ~ Li(N) - (1/2)Li(sqrt N).
So ratio = share*2h ~ 1 - (1/2)(h + a - 1) Li(sqrt N)/Li(N). Computed here independently of the page (numpy marking)."""
import numpy as np, math
from scipy.special import expi
from math import gcd
Li = lambda x: expi(math.log(x)) - expi(math.log(2))
def forms(D):
    h = a = 0
    for A in range(1, int((-D / 3) ** .5) + 1):
        for b in range(-A + 1, A + 1):
            if (b * b - D) % (4 * A): continue
            c = (b * b - D) // (4 * A)
            if c < A or (b < 0 and A == c): continue
            if gcd(gcd(A, abs(b)), c) != 1: continue
            h += 1
            if b == 0 or b == A or A == c: a += 1
    return h, a
for N in (500_000, 10_000_000):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]: s[i * i::i] = False
    pr = np.flatnonzero(s)
    obs, pred, rows = [], [], []
    for P in range(1, 101):
        rep = np.zeros(N + 1, bool)
        for b in range(0, int((N / P) ** .5) + 1):
            amax = int((N - P * b * b) ** .5); a_ = np.arange(0, amax + 1); rep[a_ * a_ + P * b * b] = True
        ok = (4 * P) % pr != 0; tot = int(ok.sum()); drawn = int(rep[pr[ok]].sum())
        h, a = forms(-4 * P)
        r_obs = drawn / tot * 2 * h; r_pred = 1 - 0.5 * (h + a - 1) * Li(N ** .5) / Li(N)
        obs.append(r_obs); pred.append(r_pred); rows.append((P, h, a, drawn, tot, r_obs, r_pred))
    obs, pred = np.array(obs), np.array(pred)
    print(f'N = {N:,}: ratio share*2h: min {obs.min():.4f} max {obs.max():.4f}; below 1 in {int((obs < 1).sum())} of 100')
    print(f'   predicted 1 - (h+a-1)Li(sqrt N)/(2Li(N)): mean deficit obs {1 - obs.mean():.5f} vs pred {1 - pred.mean():.5f}; corr(obs, pred) {np.corrcoef(obs, pred)[0, 1]:.3f}')
    d = (1 - obs) / (1 - pred); print(f'   observed/predicted deficit: median {np.median(d):.2f}, quartiles {np.quantile(d, .25):.2f}..{np.quantile(d, .75):.2f}')
    for P in (1, 5, 27, 64, 89): r = rows[P - 1]; print(f'   P={P}: h={r[1]} a={r[2]} drawn {r[3]:,}/{r[4]:,} ratio {r[5]:.4f} pred {r[6]:.4f}')
