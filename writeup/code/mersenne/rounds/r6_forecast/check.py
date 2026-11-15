#!/usr/bin/env python3
"""Independent recheck of analysis.py: no sieve, no numpy.

Walks every prime below 10^9 with gmpy2.next_prime, adds w(q) = e^g log2(a q)/q in plain floats
(a = 6 for q = 1 mod 4, a = 2 otherwise; and a = 1), and records
  - S at every exponent, pi(10^8), pi(10^9), S(2^28), S(10^9), S(142869329),
  - for every start (each exponent, and 142869329) the smallest prime q with S(q) - S(start) >= -ln(1 - level).
Then recomputes u_n, coverage, the KS statistic, log scores, the #53 forecast and the expected counts,
and compares them with results.json.
"""
import gmpy2, math, heapq, json, time, sys

EG = math.exp(0.5772156649015329)
L6, L2 = math.log2(6), 1.0
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253,
        4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011,
        24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161,
        74207281, 77232917, 82589933, 136279841]
VER, TESTED = 83880077, 142869329
LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]
N = 10**9
starts = set(EXPS) | {TESTED - 1}   # 'all exponents below T tested': T itself is not
marks = sorted([10**8, 2**28, 10**9, TESTED - 1, VER - 1])

SW = 0.0; SP = 0.0; cnt = 0
Sexp = {}; SPexp = {}; wexp = {}; wPexp = {}; Smark = {}; pimark = {}
Sstart = {}          # S(start) inclusive of primes <= start
quant = {}           # (start, level) -> prime
heap = []
mi = 0
t0 = time.time()
q = gmpy2.mpz(1)
prev_start_pending = sorted(starts)
si = 0
ei = 0; nextexp = EXPS[0]
while True:
    q = gmpy2.next_prime(q)
    qi = int(q)
    # starts that lie below qi: their S is the sum over primes <= start = current SW (before adding qi)
    while si < len(prev_start_pending) and prev_start_pending[si] < qi:
        s = prev_start_pending[si]; Sstart[s] = SW
        for l in LEVELS:
            heapq.heappush(heap, (SW - math.log1p(-l), s, l))
        si += 1
    while mi < len(marks) and marks[mi] < qi:
        Smark[marks[mi]] = SW; pimark[marks[mi]] = cnt; mi += 1
    if qi > N:
        break
    a = L6 if (qi & 3) == 1 else L2
    lq = math.log2(qi)
    w = EG * (lq + a) / qi
    wp = EG * lq / qi
    SW += w; SP += wp; cnt += 1
    while heap and SW >= heap[0][0]:
        _, s, l = heapq.heappop(heap); quant[(s, l)] = qi
    if qi == nextexp:
        Sexp[qi] = SW; SPexp[qi] = SP; wexp[qi] = w; wPexp[qi] = wp
        ei += 1; nextexp = EXPS[ei] if ei < len(EXPS) else -1
print(f"walk done: pi(1e9) = {cnt}, pi(1e8) = {pimark[10**8]}, time {time.time()-t0:.0f}s", flush=True)
assert pimark[10**8] == 5761455 and cnt == 50847534

R = json.load(open('results.json'))
bad = 0
def cmp(name, mine, theirs, tol):
    global bad
    ok = abs(mine - theirs) <= tol
    if not ok: bad += 1
    print(f"{'OK ' if ok else 'BAD'} {name}: check {mine:.6g}  analysis {theirs:.6g}")

# transitions
us = []; mids = []; c50 = c90 = 0; bitsW = bitsP = 0.0; n_ver = 0
for i, t in enumerate(R['trans']):
    p, nx = EXPS[i], EXPS[i + 1]
    assert t['p'] == p and t['q'] == nx
    dS = Sexp[nx] - Sexp[p]; u = 1 - math.exp(-dS); ulo = 1 - math.exp(-(dS - wexp[nx]))
    if abs(u - t['u']) > 1e-9 or abs(ulo - t['ulo']) > 1e-9:
        bad += 1; print('BAD u', p, u, t['u'])
    for l, Qa in zip(LEVELS, t['Q']):
        Qc = quant[(p, l)]
        tol = 0 if Qc <= 10**6 else 1000
        if abs(Qc - Qa) > tol:
            bad += 1; print('BAD quantile', p, l, Qc, Qa)
    if nx <= VER:
        n_ver += 1
        us.append(u); mids.append((u + ulo) / 2)
        c50 += (u >= 0.25 and ulo < 0.75); c90 += (u >= 0.05 and ulo < 0.95)
        bitsW += (dS - wexp[nx]) / math.log(2) - math.log2(-math.expm1(-wexp[nx]))
        dSP = SPexp[nx] - SPexp[p]
        bitsP += (dSP - wPexp[nx]) / math.log(2) - math.log2(-math.expm1(-wPexp[nx]))
print(f"u_n and all {5*len(R['trans'])} quantiles of {len(R['trans'])} transitions compared")
cmp('verified transitions', n_ver, R['coverage']['n'], 0)
cmp('coverage 50%', c50, R['coverage']['c50'], 0)
cmp('coverage 90%', c90, R['coverage']['c90'], 0)
m = sorted(mids); n = len(m)
D = max(max((k + 1) / n - x, x - k / n) for k, x in enumerate(m))
cmp('KS D (mid-PIT)', D, R['ks']['D_all'], 1e-9)
cmp('sum dS', Sexp[EXPS[n_ver]] - Sexp[2], R['ks']['sum_dS'], 1e-6)
cmp('log score W (discrete)', bitsW, R['logscore_all']['Wd'], 1e-6)
cmp('log score P a=1 (discrete)', bitsP, R['logscore_all']['Pd'], 1e-6)
hist = [0] * 10
for x in mids: hist[min(9, int(x * 10))] += 1
print('histogram', hist, 'analysis', R['ks']['hist10']); bad += hist != R['ks']['hist10']
# #53
nx = R['next53']
for l, a_, b_ in zip(LEVELS, nx['uncond'], nx['cond']):
    cmp(f'#53 uncond Q{l}', quant[(136279841, l)], a_, 1000)
    cmp(f'#53 cond   Q{l}', quant[(TESTED - 1, l)], b_, 1000)
cmp('P(none in (136279841, T])', math.exp(-(Smark[TESTED - 1] - Sexp[136279841])), nx['pnone'], 1e-5)
cmp('E (T, 2^28]', Smark[2**28] - Smark[TESTED - 1], nx['E28'], 1e-5)
cmp('E (T, 10^9]', Smark[10**9] - Smark[TESTED - 1], nx['E9'], 1e-5)
cmp('E (VER, T]', Smark[TESTED - 1] - Smark[VER - 1], nx['Ever'], 1e-5)
print(f"#53 exact prime quantiles (uncond): {[quant[(136279841, l)] for l in LEVELS]}")
print(f"#53 exact prime quantiles (cond on T): {[quant[(TESTED - 1, l)] for l in LEVELS]}")
print(f"S(10^9) = {SW:.9f}")
print('ALL OK' if bad == 0 else f'{bad} MISMATCHES')
