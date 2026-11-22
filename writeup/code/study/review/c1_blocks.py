"""Check 1: recompute block-0 statistics of windows A and C from primes listed with gmpy2.next_prime (no sieve)
and compare with results_r3.jsonl and results_twins.jsonl (exact equality required)."""
import json, gmpy2

D = '/home/bas/praime/writeup/data/study/'
BLOCKS = {'A': (10**10, 5_000_000), 'C': (10**15, 2_000_000)}


def primes_in(lo, hi):
    out = []
    p = gmpy2.next_prime(lo - 1)
    while p < hi:
        out.append(int(p))
        p = gmpy2.next_prime(p)
    return out


def stats(P):
    g = [b - a for a, b in zip(P, P[1:])]
    n = len(g)
    d = {}
    d['N.eqgap'] = sum(1 for i in range(1, n) if g[i] == g[i - 1]) / (n - 1)
    for h in (2, 4, 6, 8, 10, 12, 18, 30):
        d['S.gap%d' % h] = sum(1 for x in g if x == h) / n
    for q in (3, 4, 5, 7, 8):
        d['A.same%d' % q] = sum(1 for a, b in zip(P, P[1:]) if a % q == b % q) / (len(P) - 1)
    c3 = sum(1 for i in range(2, n) if g[i] == g[i - 1] == g[i - 2])
    d['T3.cpap4'] = c3 / n          # = mean over n-2 triples * (n-2)/n
    return d


def twin_stats(P, lo, hi):
    # P covers [lo, hi + 3)
    S = set(P)
    T = [p for p in P if lo <= p < hi and p + 2 in S and p > 5]
    h = [b - a for a, b in zip(T, T[1:])]
    d = {'n_twins': len(T), 'n_gaps': len(h)}
    for x in (6, 12, 18, 24, 30, 36, 42, 48, 60, 90, 120, 150, 180, 210):
        d['W1.gap%d' % x] = sum(1 for y in h if y == x) / len(h)
    sm = [y for y in h if y <= 300]
    d['W2.le300'] = len(sm) / len(h)
    for q in (5, 7):
        d['W3.same%d' % q] = sum(1 for y in sm if y % q == 0) / len(sm)
    return d


rec3 = {}
for l in open(D + 'results_r3.jsonl'):
    r = json.loads(l)
    if r['block'] == 0:
        rec3[r['window']] = r
rect = {}
for l in open(D + 'results_twins.jsonl'):
    r = json.loads(l)
    if r['block'] == 0:
        rect[r['window']] = r

for w, (L, W) in BLOCKS.items():
    P = primes_in(L, L + W + 3)
    Pin = [p for p in P if p < L + W]
    print(w, 'primes', len(Pin), 'stored members', rec3[w]['members'])
    mine = stats(Pin)
    bad = 0
    for k, v in mine.items():
        s = rec3[w]['stats'][k]
        ok = abs(v - s) <= 1e-15 * max(1, abs(s))
        bad += not ok
        print('  %-10s mine %.15g  stored %.15g  %s' % (k, v, s, 'OK' if ok else 'DIFF'))
    tw = twin_stats(P, L, L + W)
    for k, v in tw.items():
        s = rect[w]['stats'][k]
        ok = abs(v - s) <= 1e-15 * max(1, abs(s))
        bad += not ok
        print('  %-10s mine %.15g  stored %.15g  %s' % (k, v, s, 'OK' if ok else 'DIFF'))
    print(w, 'mismatches', bad)
