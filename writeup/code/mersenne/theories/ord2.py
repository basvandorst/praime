"""Theories 6, 7, 8, 10.
T7: the number of primes q <= x whose order of 2 is a prime l (the prime factors of 2^l - 1, l prime). Model: 2 behaves like a
random element of the group of index 2 when q = +-1 mod 8 and like a non-residue otherwise, so
P(ord = l) = 2(l - 1)/(q - 1) for q = +-1 mod 8 and every prime l | q - 1, l odd, and 0 for q = +-3 mod 8.
T10: q = 1 mod 3 exactly when 3 | k, k = (q - 1)/(2l) (proved for l != 3), so Mersenne factors lean to 2 mod 3; the model
predicts the split. T8: p = 1 mod 4, q = 6p + 1 prime: q | 2^p - 1 exactly when q = x^2 + 27y^2 (proved from Gauss's cubic
criterion), checked directly with x, y. T6: twin race, from s1e9.npz."""
import sys, math, numpy as np
X = int(float(sys.argv[1]))
spf = np.zeros(X + 1, np.int32)
for i in range(2, int(X ** 0.5) + 1):
    if spf[i] == 0:
        sl = spf[i * i::i]; sl[sl == 0] = i
isp = spf == 0; isp[:2] = False
P = np.nonzero(isp)[0]
def pf(n):
    out = []
    while n > 1:
        r = int(spf[n]) or n
        out.append(r)
        while n % r == 0: n //= r
    return out
cps = [10 ** e for e in range(3, 9) if 10 ** e <= X]
obs = {c: 0 for c in cps}; exp = {c: 0.0 for c in cps}; exp4 = {c: 0.0 for c in cps}
o3 = {1: 0, 2: 0}; e3 = {1: 0.0, 2: 0.0}
for q in P[2:]:
    q = int(q)
    fs = [l for l in pf(q - 1) if l > 2]
    hit = any(pow(2, l, q) == 1 for l in fs)
    good = q % 8 in (1, 7)
    e = sum(2 * (l - 1) / (q - 1) for l in fs) if good else 0.0
    if hit:
        o3[q % 3] += 1
    e3[q % 3] += e
    for c in cps:
        if q <= c:
            obs[c] += hit; exp[c] += e
print(f'== T7: primes q <= x with prime order of 2 (x up to {X:,}) ==')
print('       x   observed   model (2 random in the residues)   ratio')
for c in cps: print(f'{c:>8,} {obs[c]:>10,} {exp[c]:>16,.1f} {obs[c] / exp[c]:>22.4f}')
print(f'\n== T10: their residue mod 3 (q >= 5, up to {X:,}) ==')
tot = o3[1] + o3[2]; te = e3[1] + e3[2]
print(f'q = 1 mod 3: observed {o3[1]:,} ({o3[1] / tot:.4f}), model {e3[1]:,.1f} ({e3[1] / te:.4f}); q = 2 mod 3: observed {o3[2]:,}, model {e3[2]:,.1f}')
print('check of the rule q = 1 mod 3 <=> 3 | k for the order l > 3: ', end='')
bad = 0; n = 0
for q in P[2:min(len(P), 200000)]:
    q = int(q)
    for l in pf(q - 1):
        if l > 3 and pow(2, l, q) == 1:
            n += 1; k = (q - 1) // (2 * l); bad += ((q % 3 == 1) != (k % 3 == 0))
print(f'{n} factors, {bad} exceptions')
print(f'\n== T8: p = 1 mod 4, q = 6p + 1 prime, q | 2^p - 1  <=>  q = x^2 + 27 y^2 ==')
lim = (X - 1) // 6
ps = P[(P % 4 == 1) & (P <= lim)]
n = 0; both = 0; div_only = 0; form_only = 0; ndiv = 0
for p in ps:
    p = int(p); q = 6 * p + 1
    if not isp[q]: continue
    n += 1
    d = pow(2, p, q) == 1
    y = np.arange(0, int(math.isqrt(q // 27)) + 1)
    r = q - 27 * y * y; rt = np.floor(np.sqrt(r)).astype(np.int64)
    f = bool(np.any(rt * rt == r) or np.any((rt + 1) ** 2 == r))
    both += d and f; div_only += d and not f; form_only += f and not d; ndiv += d
print(f'p = 1 mod 4 up to {lim:,} with 6p + 1 prime: {n:,}; 6p + 1 divides 2^p - 1 for {ndiv:,} ({ndiv / n:.4f}, Chebotarev predicts 1/3); exceptions: divides but no form {div_only}, form but no divide {form_only}')
try:
    a = np.load('s1e9.npz')
    bu = np.cumsum(a['bu']); bd = np.cumsum(a['bd'])
    xs = 24000 * np.arange(1, len(bu) + 1)
    sp = P[(P >= 5) & (P.astype(np.float64) ** 2 <= 1.0e9 + 50)]
    import sympy
    sq2 = np.array([sympy.isprime(int(p) * int(p) - 2) for p in sp])
    Nx = np.searchsorted(sp.astype(np.float64) ** 2, xs, side='right')
    cn = np.concatenate([[0], np.cumsum(sq2)])[Nx]
    D = bd - bu
    print(f'\n== T6: twin race to 10^9: down both (12i+5, 12i+7) minus up both (12i-1, 12i+1) ==')
    print(f'twin pairs: up {bu[-1]:,}, down {bd[-1]:,}, D = {D[-1]:+,}; predicted bias +N/2 with N = #(p <= sqrt(x): p^2 - 2 prime) = {cn[-1]}, so {cn[-1] / 2:+.0f}')
    for c in (10 ** 6, 10 ** 7, 10 ** 8, 10 ** 9):
        i = np.searchsorted(xs, c) - 1
        print(f'  x {c:>13,}: D {D[i]:>+6}, predicted {cn[i] / 2:>+7.1f}, sd of a fair race sqrt(up + down) = {math.sqrt(bu[i] + bd[i]):.0f}')
    lx = np.log(xs); w = np.gradient(lx)
    sel = xs >= 1e5
    print(f'share of log-scale x in [1e5, 1e9] with D > 0: {np.sum(w[sel] * (D[sel] > 0)) / np.sum(w[sel]):.3f}; log-mean of D: {np.sum(w[sel] * D[sel]) / np.sum(w[sel]):+.1f}; log-mean of the predicted bias: {np.sum(w[sel] * cn[sel] / 2) / np.sum(w[sel]):+.1f}')
    print(f'log-mean of D/sqrt(count): {np.sum(w[sel] * D[sel] / np.sqrt(bu[sel] + bd[sel])) / np.sum(w[sel]):+.3f}; of the prediction: {np.sum(w[sel] * cn[sel] / 2 / np.sqrt(bu[sel] + bd[sel])) / np.sum(w[sel]):+.3f}')
except FileNotFoundError:
    pass
