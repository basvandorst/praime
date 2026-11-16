"""Theories 7 and 10 in closed form. A prime q with prime order l of 2 is q = 2kl + 1 with 2 a (2k)-th power residue.
For fixed k the primes l <= y with 2kl + 1 prime number about 2 C2 g(k) Li2-type(y) (Hardy-Littlewood for the pair l, 2kl + 1;
g(k) = prod over odd primes r | k of (r - 1)/(r - 2)); 2kl + 1 = +-1 mod 8 holds for every l when 4 | k, for half of the l when
k is odd, never when k = 2 mod 4 (m(k) = 1, 1/2, 0); then the order is l with chance 1/k. So
  N(x) ~ C2 * sum_k m(k) g(k) / k * 2 * int_2^{x/(2k)} dt / (ln t ln(2kt))   ~   A x / ln^2 x,   A = C2 sum_k m(k) g(k) / k^2,
and q = 1 mod 3 exactly when 3 | k, so the share of such factors tends to sum_{3|k} / sum_k (with the same weights)."""
import numpy as np, math
C2 = 0.6601618158468696
K = 2_000_000
spf = np.zeros(K + 1, np.int64)
for i in range(2, K + 1):
    if spf[i] == 0:
        sl = spf[i::i]; sl[sl == 0] = i
g = np.ones(K + 1)
for k in range(2, K + 1):
    n = k
    while n > 1:
        r = int(spf[n])
        if r > 2: g[k] *= (r - 1) / (r - 2)
        while n % r == 0: n //= r
k = np.arange(K + 1)
m = np.where(k % 4 == 0, 1.0, np.where(k % 2 == 1, 0.5, 0.0)); m[0] = 0
w = m * g
A = C2 * np.sum(w[1:] / k[1:].astype(float) ** 2)
share = np.sum(w[3::3] / k[3::3].astype(float) ** 2) / np.sum(w[1:] / k[1:].astype(float) ** 2)
print(f'A = C2 * sum m(k) g(k) / k^2 = {A:.6f};  limit share of factors = 1 mod 3: {share:.4f}')
def N(x):
    tot = 0.0; tot3 = 0.0
    for kk in range(1, int(x // 6) + 1):
        if w[kk] == 0: continue
        y = (x - 1) / (2 * kk)
        if y < 3: break
        t = np.exp(np.linspace(math.log(3), math.log(y), 4000))
        f = 1 / (np.log(t) * np.log(2 * kk * t))
        I = np.sum((f[1:] + f[:-1]) / 2 * np.diff(t))
        c = 2 * C2 * w[kk] / kk * I
        tot += c
        if kk % 3 == 0: tot3 += c
    return tot, tot3
obs = {10 ** 3: 21, 10 ** 4: 105, 10 ** 5: 585, 10 ** 6: 3845, 10 ** 7: 26560}
obs3 = {10 ** 6: 852, 10 ** 7: 5727}
for x in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    n, n3 = N(x)
    extra = f'; share 1 mod 3: formula {n3 / n:.4f}, found {obs3[x] / obs[x]:.4f}' if x in obs3 else ''
    print(f'x = {x:>10,}: formula {n:>9.1f}, found {obs[x]:>6}, ratio {obs[x] / n:.4f}; A x / ln^2 x = {A * x / math.log(x) ** 2:.1f}{extra}')
