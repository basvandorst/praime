"""Primes in short intervals (baskets of width h), windows n in [N/2, N].
(1) weighted: variance of psi(n + h) - psi(n) - h  against Montgomery-Soundararajan (2004, Theorem 3, conditional on a
    uniform Hardy-Littlewood conjecture, for log N <= h <= N^(1/2)):  h (ln(n/h) + B), B = 1 - gamma - ln 2 pi,
    averaged over the window; and against Cramer's model: h ln n averaged (density 1/ln n, weights ln n).
(2) unweighted: variance of pi(n + h) - pi(n) - (li(n + h) - li(n)) divided by the mean count; heuristic prediction
    1 - (ln h - B) / ln n (divide (1) by ln^2 n), averaged; Poisson gives 1; and the same for fakes (true primes up to
    max(W, 100), then numbers without prime factors in W, chance prod p/(p-1) / ln n), detrended by their own expectation.
Prediction written before the run: (1) and (2) match the primes for h up to sqrt(N); the fakes stay near 1 (W = {2})
and come down as W grows.  The fakes are shown twice: detrended by their own exact expectation (which removes the
sieve's own regular pattern, so the ratio is 1 - c_W / ln n for every h) and detrended by li like the primes (the fair
comparison, added after a review).   usage: baskets.py N runs out.txt"""
import sys, math
import numpy as np
N, RUNS = int(sys.argv[1]), int(sys.argv[2]); out = open(sys.argv[3], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
HS = [h for h in (10, 30, 100, 300, 1000, 3000, 10000) if h <= N // 20]
M = N + max(HS) + 1
sv = np.ones(M + 1, bool); sv[:2] = False
for i in range(2, int(M ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
lam = np.zeros(M + 1)
for p in P:
    lp = math.log(p); q = int(p)
    while q <= M: lam[q] += lp; q *= int(p)
psi = np.cumsum(lam)
B = 1 - 0.5772156649015329 - math.log(2 * math.pi)
n = np.arange(N // 2, N + 1); ln = np.log(n)
t = np.arange(M + 1, dtype=float); dens = np.zeros(M + 1); dens[2:] = 1 / np.log(t[2:]); li = np.cumsum(dens)
def ratio(ind, expc, h):
    c = np.cumsum(ind)
    d = (c[n + h] - c[n]) - (expc[n + h] - expc[n]); return float(d.var() / (expc[n + h] - expc[n]).mean())
pr('N = %d, window [N/2, N], B = %.4f, sqrt(N) = %.0f' % (N, B, math.sqrt(N)))
pr('     h | weighted: var     MS     Cramer  var/MS var/Cramer | unweighted var/mean: primes  heuristic | fakes, detrended by their own expectation: W=2 W<=13 W<=100 | fakes detrended by li like the primes: W=2 W<=13 W<=100')
pr('runs per fake model: %d' % RUNS)
fake = {}
rng = np.random.default_rng(16)
for Bw in (2, 13, 100):
    W = P[P <= Bw]; ok = np.ones(M + 1, bool); ok[:max(Bw, 100) + 1] = False
    for p in W: ok[::p] = False
    c = float(np.prod(W / (W - 1.0))); pb = np.zeros(M + 1); pb[ok] = np.minimum(1, c / np.log(t[ok])); pb[P[P <= max(Bw, 100)]] = 1
    fake[Bw] = (pb, [(rng.random(M + 1) < pb).astype(float) for _ in range(RUNS)])
for h in HS:
    d = psi[n + h] - psi[n] - h; v = float(d.var()); ms = float(np.mean(h * (np.log(n / h) + B))); cr = float(np.mean(h * ln))
    rp = ratio(sv.astype(float), li, h); heur = float(np.mean(1 - (math.log(h) - B) / ln))
    fr = [np.mean([ratio(f, np.cumsum(pb), h) for f in fs]) for pb, fs in (fake[2], fake[13], fake[100])]
    fl = [np.mean([ratio(f, li, h) for f in fs]) for pb, fs in (fake[2], fake[13], fake[100])]
    pr('%6d | %13.1f %8.1f %8.1f %6.3f %8.3f | %24.3f %9.3f | %9.3f %7.3f %7.3f | %9.3f %7.3f %7.3f' % (h, v, ms, cr, v / ms, v / cr, rp, heur, *fr, *fl))
