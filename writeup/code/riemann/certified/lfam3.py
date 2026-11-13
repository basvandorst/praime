"""Zeros of all Dirichlet L-functions mod a prime q on 0 < t <= T, for any height (lfam.py with an adaptive
Euler-Maclaurin length N(t) so that (|s| + 2M)/(2 pi N) <= 0.3), the full grid of Z kept in memory, the argument
principle per character, and a repair pass: a minimum of |Z| without a sign change marks a close pair inside one grid
step; it is resolved by direct evaluation. Usage: lfam3.py q T h out.npz"""
import sys, math, time
import numpy as np
from scipy.special import loggamma, bernoulli
from sympy import isprime, primitive_root
q = int(sys.argv[1]); T = float(sys.argv[2]); h = float(sys.argv[3]); out = sys.argv[4]
assert isprime(q)
g = int(primitive_root(q)); pw = np.empty(q - 1, dtype=np.int64); v = 1
for j in range(q - 1): pw[j] = v; v = v * g % q
alpha = pw / q; M_ = 9
Nmax = max(30, int(math.ceil((T + 2 * M_ + 2) / (0.6 * math.pi))) + 2)
LOGALL = np.log(np.arange(Nmax)[None, :] + alpha[:, None])
B = bernoulli(2 * M_); coef = [float(B[2 * j]) / math.factorial(2 * j) for j in range(1, M_ + 1)]
def hurwitz_all(s):
    N = max(30, int(math.ceil((abs(s) + 2 * M_) / (0.6 * math.pi))))
    main = np.exp(-s * LOGALL[:, :N]).sum(axis=1)
    Na = N + alpha; Nas = np.exp(-s * np.log(Na))
    tail = Nas * Na / (s - 1) + 0.5 * Nas; rising = s; powr = Nas / Na
    for j in range(1, M_ + 1):
        tail = tail + coef[j - 1] * rising * powr; rising = rising * (s + 2 * j - 1) * (s + 2 * j); powr = powr / (Na * Na)
    return main + tail
def L_all(s): return np.exp(-s * math.log(q)) * (q - 1) * np.fft.ifft(hurwitz_all(s))
k = np.arange(q - 1); kappa = k % 2; kbar = (-k) % (q - 1)
tau = (q - 1) * np.fft.ifft(np.exp(2j * math.pi * pw / q)); eps = tau / ((1j) ** kappa * math.sqrt(q))
rot = np.exp(-0.5j * np.angle(eps))
def theta(t, kap): return 0.5 * t * math.log(q / math.pi) + np.imag(loggamma((0.5 + kap + 1j * t) / 2))
def Zrow(t):
    W = rot * np.exp(1j * np.where(kappa == 0, theta(t, 0), theta(t, 1))) * L_all(0.5 + 1j * t)
    return W.real, float(np.max(np.abs(W.imag[1:]) / (np.abs(W[1:]) + 1e-3)))
t0 = time.time(); pad = 8
ts = (np.arange(int(math.ceil(T / h)) + 2 * pad + 1) - pad) * h
Z = np.empty((len(ts), q - 1)); imrat = 0.0
for i, t in enumerate(ts):
    Z[i], r = Zrow(t); imrat = max(imrat, r)
print(f'grid: {len(ts)} rows up to t = {ts[-1]:.2f}, N up to {Nmax}, {time.time() - t0:.0f} s; max |Im W|/(|W| + 1e-3) = {imrat:.2e}', flush=True)
U = np.arange(-3, 5, dtype=float); Vinv = np.linalg.inv(np.vander(U, 8, increasing=True))
U5 = np.arange(-2, 4, dtype=float); V5 = np.linalg.inv(np.vander(U5, 6, increasing=True))
def newton(C, u):
    for _ in range(40):
        p = np.zeros_like(u); dp = np.zeros_like(u)
        for c_ in C.T[::-1]: dp = dp * u + p; p = p * u + c_
        u = np.clip(u - p / dp, 0, 1)
    return u
rows, cols = np.nonzero(Z[:-1] * Z[1:] < 0)
sel = (rows >= 3) & (rows + 4 < len(ts)) & (cols >= 1); rows, cols = rows[sel], cols[sel]
vals = np.stack([Z[rows + d, cols] for d in range(-3, 5)], axis=1)
u = newton(vals @ Vinv.T, Z[rows, cols] / (Z[rows, cols] - Z[rows + 1, cols]))
u5 = newton(vals[:, 1:7] @ V5.T, u.copy())
gam = ts[rows] + u * h; interr = np.abs(u5 - u) * h
ok = (gam > 0) & (gam <= T); gk, gg, ge = cols[ok], gam[ok], interr[ok]
print(f'{len(gg)} zeros on (0, T]; interpolation difference deg 7 vs 5: max {interr.max():.2e}, median {np.median(interr):.2e}', flush=True)
# argument principle at T
sig = np.concatenate([np.linspace(2.0, 0.6, 141), np.linspace(0.6, 0.5, 201)[1:]])
prev = np.angle(L_all(sig[0] + 1j * T)); A = prev.copy(); maxstep = 0.0
for s_ in sig[1:]:
    cur = np.angle(L_all(s_ + 1j * T)); d = (cur - prev + math.pi) % (2 * math.pi) - math.pi
    maxstep = max(maxstep, float(np.abs(d[1:]).max())); A += d; prev = cur
NAP = np.rint((2 * np.where(kappa == 0, theta(T, 0), theta(T, 1)) + A + A[kbar]) / math.pi).astype(int)
zero_at_0 = np.abs(Z[pad]) < 1e-9
def counts(gk_):
    f = np.bincount(gk_, minlength=q - 1); return f + f[kbar] + zero_at_0
mism = np.nonzero((NAP != counts(gk)) & (k >= 1))[0]
print(f'argument principle: largest step {maxstep:.3f} rad; characters whose count differs: {len(mism)}', flush=True)
# repair: dips of |Z| without a sign change, for the mismatched characters and their conjugates
add_k, add_g = [], []
jj = np.arange(q - 1)
for kk in sorted(set(mism.tolist()) | set(kbar[mism].tolist())):
    z = Z[:, kk]; ch = np.exp(2j * math.pi * ((kk * jj) % (q - 1)) / (q - 1))
    def Zk(t):
        s = 0.5 + 1j * t
        Lk = np.exp(-s * math.log(q)) * (ch @ hurwitz_all(s))
        return float(np.real(rot[kk] * np.exp(1j * theta(t, kappa[kk])) * Lk))
    for i in range(1, len(ts) - 1):
        if z[i - 1] * z[i] > 0 and z[i] * z[i + 1] > 0 and abs(z[i]) < abs(z[i - 1]) and abs(z[i]) < abs(z[i + 1]) and 0 < ts[i] <= T + h:
            den = z[i - 1] - 2 * z[i] + z[i + 1]; tv = ts[i] + 0.5 * h * (z[i - 1] - z[i + 1]) / den
            zv = Zk(tv)
            if zv * z[i] < 0:
                for lo, hi in ((min(ts[i - 1], tv), max(ts[i - 1], tv)) if tv < ts[i] else (ts[i], tv), (tv, ts[i + 1]) if tv > ts[i] else (tv, ts[i])):
                    zlo = Zk(lo)
                    if zlo * Zk(hi) > 0: continue
                    for _ in range(50):
                        mid = 0.5 * (lo + hi); zm = Zk(mid)
                        if zm * zlo > 0: lo, zlo = mid, zm
                        else: hi = mid
                    gm = 0.5 * (lo + hi)
                    if 0 < gm <= T: add_k.append(kk); add_g.append(gm)
gk = np.concatenate([gk, np.array(add_k, dtype=gk.dtype)]); gg = np.concatenate([gg, np.array(add_g)])
ge = np.concatenate([ge, np.full(len(add_g), 1e-10)])
cnt = counts(gk); mism2 = np.nonzero((NAP != cnt) & (k >= 1))[0]
print(f'repair: added {len(add_g)} zeros; characters whose count still differs: {len(mism2)} {mism2[:10].tolist()}; {time.time() - t0:.0f} s', flush=True)
np.savez(out, q=q, g=g, T=T, h=h, k=gk, gamma=gg, interr=ge, count=cnt, NAP=NAP, pw=pw, mism=mism2)
