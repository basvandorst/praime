"""P-1 on M_p = 2^p - 1 for the trial-factoring survivors, run for real with gmpy2.
  pm1.py X K lo hi B1levels S2B1 B2checkpoints out.json
Survivors: primes lo <= p < hi with no factor q = 2kp + 1, k <= K (kmin file kmin_<X>_<K'>.npy, K' >= K).
Stage 1: x = 3^E mod M_p with E = 2p * prod_{l <= B1} l^floor(log_l B1), run through the nested B1 levels
(x is raised by E_B / E_Bprev); gcd(x - 1, M_p) after each level. Cost = bitlen(E_B) squarings.
Stage 2 (standard, prime by prime): from x_S = x at B1 = S2B1, for every prime r in (S2B1, B2]: y = x_S^r (updated
by the precomputed x_S^gap), acc *= (y - 1); gcd(acc, M_p) at each checkpoint. Cost = 2 multiplications per prime.
A factor g is recorded only if 1 < g < M_p (then M_p is proved composite); g = M_p is inconclusive (tiny p only).
Output per p: first B1 level that finds a factor, first B2 checkpoint that finds one (if stage 1 at S2B1 did not), g."""
import sys, json, time, math, numpy as np, gmpy2

X, K, lo, hi = map(int, sys.argv[1:5])
LEV = [int(v) for v in sys.argv[5].split(',')]
S2B1 = int(sys.argv[6]); CP = [int(v) for v in sys.argv[7].split(',')]
OUT = sys.argv[8]
assert S2B1 in LEV
files = {10000: [10000000, 1000000], 100000: [1000000, 10000]}[X]
for Kf in files:
    try: P, kmin = np.load(f'kmin_{X}_{Kf}.npy'); break
    except FileNotFoundError: P, kmin = np.load(f'../kmin_{X}_{Kf}.npy')
assert Kf >= K
surv = [int(p) for p, k in zip(P, kmin) if lo <= p < hi and (k == 0 or k > K)]
STRIDE = int(sys.argv[9]) if len(sys.argv) > 9 else 1      # optional: every STRIDE-th survivor only (a sample)
surv = surv[::STRIDE]

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]

PR = primes_upto(max(CP))
def Epow(B):                       # prod of l^floor(log_l B), l <= B
    e = gmpy2.mpz(1)
    for l in PR:
        if l > B: break
        a = l
        while a * l <= B: a *= l
        e *= a
    return e
EB = {B: Epow(B) for B in LEV}
S1COST = {B: int(gmpy2.bit_length(EB[B])) for B in LEV}            # + bitlen(2p), added per p
s2pr = [r for r in PR if S2B1 < r <= max(CP)]
S2COST = {c: 2 * sum(1 for r in s2pr if r <= c) for c in CP}

res = {}
t0 = time.time()
for n, p in enumerate(surv):
    M = (gmpy2.mpz(1) << p) - 1
    def mulmod(a, b):
        t = a * b; t = (t & M) + (t >> p)
        return t - M if t >= M else t
    x = gmpy2.powmod(3, 2 * p, M); prevE = gmpy2.mpz(1)
    r = {'s1': None, 's2': None, 'g': None, 'xS': None}
    for B in LEV:
        x = gmpy2.powmod(x, EB[B] // prevE, M); prevE = EB[B]
        g = gmpy2.gcd(x - 1, M)
        if 1 < g < M:
            r['s1'] = B; r['g'] = str(g); break
        if B == S2B1: xS = x
    if r['s1'] is None or r['s1'] > S2B1:
        # stage 2 from x_S (stage 1 at S2B1 found nothing)
        gaps = {}
        y = gmpy2.powmod(xS, s2pr[0], M); acc = (y - 1) % M; last = s2pr[0]
        ci = 0
        for rr in s2pr[1:] + [None]:
            while ci < len(CP) and (rr is None or rr > CP[ci]):
                g = gmpy2.gcd(acc, M)
                if 1 < g < M:
                    r['s2'] = CP[ci]
                    if r['g'] is None: r['g'] = str(g)
                    break
                ci += 1
            if r['s2'] is not None or rr is None: break
            d = rr - last
            if d not in gaps: gaps[d] = gmpy2.powmod(xS, d, M)
            y = mulmod(y, gaps[d]); acc = mulmod(acc, y - 1); last = rr
    del r['xS']
    res[p] = r
    if n % 100 == 0: print(f'{n}/{len(surv)} p={p} {time.time() - t0:.0f} s', flush=True)
json.dump({'stride': STRIDE, 'X': X, 'K': K, 'lo': lo, 'hi': hi, 'levels': LEV, 'S2B1': S2B1, 'checkpoints': CP,
           'S1COST': S1COST, 'S2COST': S2COST, 'res': res, 'secs': time.time() - t0}, open(OUT, 'w'))
f1 = {B: sum(1 for v in res.values() if v['s1'] is not None and v['s1'] <= B) for B in LEV}
f2 = {c: sum(1 for v in res.values() if (v['s1'] is not None and v['s1'] <= S2B1) or (v['s2'] is not None and v['s2'] <= c)) for c in CP}
print(f'{len(surv)} survivors of TF k <= {K} in [{lo},{hi}); stage 1 removes (cumulative) {f1}; with stage 2 from {S2B1}: {f2}; {time.time() - t0:.0f} s')
