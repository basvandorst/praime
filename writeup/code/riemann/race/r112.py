"""Round 112 (preregister/round112.txt): the race mod 4 from the zeros of L(s, chi_-4) only: first sign change of
E(x) = 1 + 2 sum Re(x^(i g)/rho) and the prime-weighted share of x <= 1.4e11 with E < 0. usage: r112.py"""
import os, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
g = np.array([float(l) for l in open(os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'D-4.txt')) if l.strip()])
rho = 0.5 + 1j * g; w = 2 / rho
u0, u1 = math.log(100), math.log(1.4e11)
du = 2 * math.pi / g.max() / 20
u = np.arange(u0, u1, du)
E = np.ones_like(u)
for i0 in range(0, len(u), 20000):
    uu = u[i0:i0 + 20000]
    E[i0:i0 + 20000] += np.real(np.exp(1j * np.outer(uu, g)) @ w)
neg = E < 0
first = float(np.exp(u[np.argmax(neg)])) if neg.any() else None
wt = np.exp(u) / u * du
share = float(np.sum(wt[neg]) / np.sum(wt))
# the sign-change intervals
idx = np.nonzero(np.diff(neg.astype(int)))[0]
iv = [(float(np.exp(u[a + 1])), float(np.exp(u[b + 1]))) for a, b in zip(idx[::2], idx[1::2])][:12]
out = {'zeros': int(len(g)), 'T': float(g.max()), 'grid_step_u': du, 'first_E_negative_x': first, 'share_primes_E_negative': share, 'intervals': iv}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(H, 'r112.json'), 'w'), indent=1)
