import sys, time
import numpy as np
from math import log, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, segments, zeros, spectrum, main_term, primes_upto

t0 = time.time()
X = 10**9
andrica_max = (0.0, 0)
leg_min = (10**9, 0)
brocard_min = (10**9, 0)
prev = None
sq_counts = {}
for lo, P in segments(X, 10**8):
    idx = np.flatnonzero(P[:min(10**8, X - lo)]) + lo
    if prev is not None:
        idx = np.concatenate([[prev], idx])
    g = np.diff(idx)
    a = np.sqrt(idx[1:].astype(np.float64)) - np.sqrt(idx[:-1].astype(np.float64))
    i = a.argmax()
    if a[i] > andrica_max[0]:
        andrica_max = (a[i], int(idx[i]))
    prev = idx[-1]
    n_lo, n_hi = int(np.ceil(sqrt(lo))), int(sqrt(lo + 10**8))
    ns = np.arange(max(n_lo, 1), n_hi)
    cnt = np.searchsorted(idx, (ns + 1) ** 2) - np.searchsorted(idx, ns ** 2)
    ok = (ns ** 2 >= max(lo, 2)) & ((ns + 1) ** 2 <= idx[-1])
    if ok.any() and cnt[ok].min() < leg_min[0]:
        j = np.argmin(np.where(ok, cnt, 10**9)); leg_min = (int(cnt[j]), int(ns[j]))
    print(f"  segment {lo:.0e}  [{time.time()-t0:.0f} s]", flush=True)
print(f"E-C15 Andrica to 1e9: max sqrt(p_(n+1)) - sqrt(p_n) = {andrica_max[0]:.4f} at p = {andrica_max[1]} (the conjecture says < 1; the maximum is at p = 7, 0.6708)")
print(f"E-C16 Legendre to 1e9: least number of primes in (n^2, (n+1)^2) over 1 <= n < 31622: {leg_min[0]} at n = {leg_min[1]} (expected minimum 2 at n = 5 and 9, FACTS 3020)")

print("E-F7 return map of normalised gaps near 1e8: corr(g_n, g_(n+1)) and corr(g_n, g_(n+2)) against the exponential model (0) and the sieve model (FINDINGS: -0.029):")
isp = sieve(10**8)
pr = np.flatnonzero(isp[10**7:10**8]) + 10**7
g = np.diff(pr).astype(np.float64)
print(f"  corr lag 1: {np.corrcoef(g[:-1], g[1:])[0,1]:+.4f}, lag 2: {np.corrcoef(g[:-2], g[2:])[0,1]:+.4f}, lag 3: {np.corrcoef(g[:-3], g[3:])[0,1]:+.4f}; P(g_(n+1) = g_n) = {(g[1:] == g[:-1]).mean():.4f} against {(np.random.default_rng(0).permutation(g)[1:] == g[:-1]).mean():.4f} for shuffled gaps")

print("E-Y1 Gilbreath: iterated absolute differences of the first 100000 primes, leading entries of rows 1..5000 all 1? and the row at which all entries after the first are <= 2:")
P = primes_upto(1299709)[:100000].astype(np.int64)
row = P.copy()
first_bad = None
small_row = None
for k in range(1, 5001):
    row = np.abs(np.diff(row))
    if row[0] != 1 and first_bad is None:
        first_bad = k
    if small_row is None and (row[1:] <= 2).all():
        small_row = k
print(f"  rows 1..5000: first row whose leading entry is not 1: {first_bad}; first row where every other entry is 0 or 2: {small_row} (Odlyzko's heuristic: once a row is all 0, 2 beyond the first entry, the leading 1 persists)")

print("E-R4 Linnik: the least prime P(a, q) = a mod q over all a coprime to q, for q <= 5000 (primes to 5e6): max of P(a,q)/(phi(q) ln^2 q) and whether P(a,q) < q^2 always:")
pr_all = np.flatnonzero(isp)
pr_lin = pr_all[pr_all <= 5 * 10**6]
worst = (0.0, 0, 0)
viol = 0
for q in range(3, 5001):
    r = pr_lin % q
    first = np.full(q, -1)
    seen = np.zeros(q, dtype=bool)
    order = np.argsort(r, kind='stable')
    rs = r[order]
    uniq, pos = np.unique(rs, return_index=True)
    first[uniq] = pr_lin[order][pos]
    cop = np.array([a for a in range(1, q) if np.gcd(a, q) == 1])
    Pq = first[cop].max()
    phi = len(cop)
    ratio = Pq / (phi * log(q) ** 2)
    if Pq >= q * q:
        viol += 1
    if ratio > worst[0]:
        worst = (ratio, q, int(Pq))
    if q % 2500 == 0:
        print(f"    q <= {q}: worst ratio so far {worst[0]:.3f} at q = {worst[1]} (P = {worst[2]}), violations of P < q^2: {viol}  [{time.time()-t0:.0f} s]", flush=True)
print(f"  max P(a,q)/(phi(q) ln^2 q) = {worst[0]:.3f} at q = {worst[1]} (P = {worst[2]}); cases with P(a,q) >= q^2 for 3 <= q <= 5000: {viol} (Linnik's constant L = 2 conjecture)")

print("E-A12 blind zero finding accuracy against the size of the prime set: rms distance of the recovered first 13 zeros to Odlyzko for primes to 1e5, 1e6, 1e7 (p >= 1000):")
gam = zeros(60)
true = gam[(gam > 10) & (gam < 60)]
tg = np.arange(10.0, 60.0, 0.02)
for Xb in (10**5, 10**6, 10**7):
    pb = pr_all[(pr_all >= 1000) & (pr_all <= Xb)].astype(np.float64)
    u = log(Xb / 1000)
    S = (spectrum(pb, np.log(pb), tg) - main_term(Xb, 1000, tg)).real / u
    found = []
    for i in range(1, len(tg) - 1):
        if S[i] < S[i - 1] and S[i] < S[i + 1] and S[i] < -0.5:
            a, b, c = S[i - 1], S[i], S[i + 1]
            found.append(tg[i] + 0.5 * (a - c) / (a - 2 * b + c) * 0.02)
    found = np.array(found)
    d = np.array([np.abs(found - g).min() for g in true]) if len(found) else np.full(len(true), 9.0)
    print(f"  X = {Xb:.0e}: {len(found)} dips, zeros within 0.1: {(d < 0.1).sum()}/{len(true)}, rms distance {sqrt((d[d < 0.5]**2).mean()) if (d < 0.5).any() else float('nan'):.4f}, spurious {sum(np.abs(true - f).min() > 0.1 for f in found)}")
print(f"[{time.time()-t0:.0f} s]")
