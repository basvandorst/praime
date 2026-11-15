"""Part 3: oracle-free ordering of the LL tests among primes p < 10^5, and a test of the posterior formula.
Information available without an LL test: p, the TF depth K (no factor 2kp+1, k <= K) and the P-1 outcome
(B1 = 1000, B2 = 3000 from pm1_*.json).  Candidates = survivors of TF and P-1 (p = 2 included).
Posteriors (heuristic):
  prior (Wagstaff)      pi0 = e^g log2(a p)/p, a = 2 (p = 3 mod 4) or 6 (p = 1 mod 4);
  after TF to Q = 2Kp+1  pi  = pi0 / P(no factor <= Q) = e^g log2(Q)/p  (a cancels: P(no factor <= Q) ~ log(ap)/log Q);
  after P-1 too          pi1 = pi / (1 - f), f = fraction of TF survivors that P-1 factored in the same p-decade bucket.
Orderings compared: ascending p; descending pi1; descending pi0 (prior with the a factor, ignores TF depth).
Also: does the a factor survive TF?  Removal fraction of TF from k = 10^4 to 10^6 by p mod 4 vs 1 - ln(Q1)/ln(Q2)."""
import json, math, numpy as np
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243]
E = set(MEXP); eg = math.exp(0.5772156649015329)
P, km = np.load('kmin_100000_1000000.npy'); P = [int(p) for p in P]; km = [int(k) for k in km]
pm = {}
for f in ['pm1_1e4.json', 'pm1_1e5a.json', 'pm1_1e5b.json', 'pm1_1e5c.json']:   # c is a 1-in-3 sample
    d = json.load(open(f))
    for p, r in d['res'].items():
        found = (r['s1'] is not None and r['s1'] <= 1000) or (r['s2'] is not None and r['s2'] <= 3000)
        pm[int(p)] = found

# --- does p mod 4 matter after TF?  factors with 10^4 < k <= 10^6 among survivors of k <= 10^4
print('TF removal in 10^4 < k <= 10^6 among survivors of k <= 10^4, p > 1000:')
for cls in (1, 3):
    s = [(p, k) for p, k in zip(P, km) if p > 1000 and p % 4 == cls and (k == 0 or k > 10 ** 4)]
    rem = sum(1 for p, k in s if k)
    pred = sum(1 - math.log(2e4 * p) / math.log(2e6 * p) for p, k in s)
    print(f'  p = {cls} mod 4: {len(s)} survivors, removed {rem} ({rem / len(s):.1%}), heuristic {pred:.0f} ({pred / len(s):.1%})')
s = [(p, k) for p, k in zip(P, km) if p > 1000]
for cls in (1, 3):
    n = sum(1 for p, k in s if p % 4 == cls); r = sum(1 for p, k in s if p % 4 == cls and 0 < k <= 10 ** 4)
    print(f'  (prior side) p = {cls} mod 4, p > 1000: {n} primes, factor with k <= 10^4: {r / n:.1%}')

def run(K, use_pm1):
    cand = [2] + [p for p, k in zip(P, km) if (k == 0 or k > K) and not (use_pm1 and pm.get(p, False))]
    # P-1 boost per decade bucket (computed from the cheap-test outcomes only)
    boost = {}
    if use_pm1:
        for b in range(0, 6):
            lo, hi = 10 ** b, 10 ** (b + 1)
            tfs = [p for p, k in zip(P, km) if lo <= p < hi and (k == 0 or k > K) and p in pm]
            boost[b] = 1 / (1 - sum(pm[p] for p in tfs) / len(tfs)) if tfs else 1.0
    def pi(p):
        if p == 2: return 1.0
        b = boost.get(int(math.log10(p)), 1.0)
        return min(1.0, eg * math.log2(2 * K * p + 1) / p * b)
    def pi0(p):
        a = 6 if p % 4 == 1 else 2
        return min(1.0, eg * math.log2(a * p) / p)
    orders = {'ascending': sorted(cand), 'by posterior': sorted(cand, key=lambda p: -pi(p)), 'by prior (a factor)': sorted(cand, key=lambda p: -pi0(p))}
    if use_pm1:
        cs = [p for p, k in zip(P, km) if 70000 <= p and (k == 0 or k > K)]
        smp = [p for p in cs if p in pm]; fr = sum(pm[p] for p in smp) / len(smp)
        print(f'\n  [7e4,1e5): {len(cs)} TF survivors, P-1 run on {len(smp)} (sample), factored {fr:.1%}; '
              f'estimated further removals among the {len(cs) - len(smp)} unsampled: {fr * (len(cs) - len(smp)):.0f} (counted as LL tests below)')
    print(f'\nK = {K}, P-1 {"(B1 1000, B2 3000)" if use_pm1 else "none"}: {len(cand)} LL tests needed for p < 10^5, '
          f'expected exponents sum(pi) = {sum(pi(p) for p in cand):.1f}, actual {sum(p in E for p in cand)}')
    Ns = [10, 20, 50, 100, 200, 500, 1000, 2000, len(cand)]
    print('  n tests:' + ''.join(f'{n:>7}' for n in Ns))
    for name, o in orders.items():
        act = [sum(p in E for p in o[:n]) for n in Ns]
        exp_ = [sum(pi(p) for p in o[:n]) for n in Ns]
        print(f'  {name:<20} found ' + ' '.join(f'{a:>6}' for a in act))
        print(f'  {"":<20} expect' + ' '.join(f'{x:>6.1f}' for x in exp_))
    # first-m exponents: LL tests needed
    for name, o in orders.items():
        pos = [i + 1 for i, p in enumerate(o) if p in E]
        print(f'  {name:<20} tests to find the first 10/20/25/28 exponents: {pos[9]} {pos[19]} {pos[24]} {pos[27]}')
    same = sum(1 for a, b in zip(orders['ascending'], orders['by posterior']) if a != b)
    print(f'  positions where the posterior order differs from ascending: {same}')

for K in (10 ** 4, 10 ** 6):
    run(K, False)
run(10 ** 6, True)
