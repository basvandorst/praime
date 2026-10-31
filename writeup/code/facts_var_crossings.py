# Var/mean of #primes in windows of length H near 1e9: real vs sieve-to-x^c + coin, for the Buchstab crossings c.
import numpy as np, math
from facts_maxgap_is import PR, logprod, GAMMA
X = 10**9; L = 10**8
def omega_num(u):
    h = 1e-4; us = np.arange(1, u + h, h); w = np.where(us <= 2, 1 / us, 0.0); n1 = int(round(1 / h))
    for i in range(len(us)):
        if us[i] > 2: w[i] = (us[i-1] * w[i-1] + h * 0.5 * (w[i-1-n1] + w[i-n1])) / us[i]
    return w[-1]
def sieve(y):
    s = np.ones(L, bool)
    for q in PR[PR <= y]: q = int(q); s[(-X) % q::q] = False
    return s
Hs = [100, 1000, 10**4, 10**5, 10**6]
real = sieve(math.isqrt(X + L))
def varmean(arr_counts):
    # detrend linearly, return var/mean
    t = np.arange(len(arr_counts)); a, b = np.polyfit(t, arr_counts, 1)
    return np.var(arr_counts - (a * t + b)) / arr_counts.mean()
print("H      " + "  ".join(f"{H:>8d}" for H in Hs))
print("real   " + "  ".join(f"{varmean(real[:L//H*H].reshape(-1, H).sum(1).astype(float)):8.3f}" for H in Hs), flush=True)
for c in (0.5615, 0.4027, 0.3099):
    y = X ** c; S = sieve(y)
    rho = 1.0 / (math.log(X) * math.exp(logprod(1, y)) * math.exp(GAMMA) * omega_num(1 / c))
    row = []
    for H in Hs:
        s = S[:L//H*H].reshape(-1, H).sum(1).astype(float)
        t = np.arange(len(s)); a, b = np.polyfit(t, s, 1); vs = np.var(s - (a * t + b))
        row.append((rho**2 * vs + rho * (1 - rho) * s.mean()) / (rho * s.mean()))
    print(f"c={c}  " + "  ".join(f"{r:8.3f}" for r in row) + f"   (rho={rho:.3f})", flush=True)
lx = math.log(X)
print("MS     " + "  ".join(f"{1 - (math.log(H) + 1 + GAMMA + math.log(2*math.pi)) / lx:8.3f}" for H in Hs))
