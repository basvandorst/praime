import sympy, math, numpy as np
X = 20000.0; d = 0.1
lam = np.zeros(int(1.1*X) + 5)
for p in sympy.primerange(2, len(lam)):
    v = p
    while v < len(lam): lam[v] += math.log(p); v *= p
psi = np.cumsum(lam)
N = 4_000_000; x = 1 + (np.arange(N) + 0.5) * (X - 1) / N
f = (psi[np.floor(1.1 * x).astype(int)] - psi[np.floor(x).astype(int)] - d * x) ** 2
print('brute R', f.sum() * (X - 1) / N / (0.5 * d * X * X))
