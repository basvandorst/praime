"""Streaming version of lfam.py for large q: zeros of all Dirichlet L-functions mod a prime q on 0 < t <= T.

Same mathematics as lfam.py (Euler-Maclaurin Hurwitz values for every a/q, one FFT in discrete-log order, real
Z_k(t) = Re(eps_k^(-1/2) e^{i theta_k(t)} L(1/2+it, chi_k)), sign changes on a grid of step h refined by a
degree-7 interpolant, argument-principle count per character). Only 8 grid rows are kept in memory.
Usage: lfam2.py q T h out.npz
"""
import sys, math, time
from collections import deque
import numpy as np
from scipy.special import loggamma, bernoulli
from sympy import isprime, primitive_root

q = int(sys.argv[1]); T = float(sys.argv[2]); h = float(sys.argv[3]); out = sys.argv[4]
assert isprime(q)
g = int(primitive_root(q))
pw = np.empty(q - 1, dtype=np.int64); v = 1
for j in range(q - 1):
    pw[j] = v; v = v * g % q
alpha = pw / q
NT_, M_ = 30, 9
LOG = np.log(np.arange(NT_)[None, :] + alpha[:, None])
Na = NT_ + alpha; LNa = np.log(Na)
B = bernoulli(2 * M_)
coef = [float(B[2 * j]) / math.factorial(2 * j) for j in range(1, M_ + 1)]

def hurwitz_all(s):
    main = np.exp(-s * LOG).sum(axis=1)
    Nas = np.exp(-s * LNa)
    tail = Nas * Na / (s - 1) + 0.5 * Nas
    rising = s; powr = Nas / Na
    for j in range(1, M_ + 1):
        tail = tail + coef[j - 1] * rising * powr
        rising = rising * (s + 2 * j - 1) * (s + 2 * j)
        powr = powr / (Na * Na)
    return main + tail

def L_all(s):
    return np.exp(-s * math.log(q)) * (q - 1) * np.fft.ifft(hurwitz_all(s))

k = np.arange(q - 1); kappa = k % 2
tau = (q - 1) * np.fft.ifft(np.exp(2j * math.pi * pw / q))
eps = tau / ((1j) ** kappa * math.sqrt(q))
print('|eps| range (k >= 1):', np.abs(eps[1:]).min(), np.abs(eps[1:]).max(), flush=True)
rot = np.exp(-0.5j * np.angle(eps))
def theta(t, kap):
    return 0.5 * t * math.log(q / math.pi) + np.imag(loggamma((0.5 + kap + 1j * t) / 2))

U7 = np.arange(-3, 5, dtype=float); V7 = np.linalg.inv(np.vander(U7, 8, increasing=True))
U5 = np.arange(-2, 4, dtype=float); V5 = np.linalg.inv(np.vander(U5, 6, increasing=True))
def newton(C, u):
    for _ in range(30):
        p = np.zeros_like(u); dp = np.zeros_like(u)
        for c_ in C.T[::-1]:
            dp = dp * u + p; p = p * u + c_
        u = np.clip(u - p / dp, 0, 1)
    return u

pad = 8
ts = (np.arange(int(math.ceil(T / h)) + 2 * pad + 1) - pad) * h
buf = deque(maxlen=8); rowid = deque(maxlen=8)
KZ, GZ, ERR = [], [], []
imrat = 0.0; z0row = None
t0 = time.time()
for i, t in enumerate(ts):
    L = L_all(0.5 + 1j * t)
    th = np.where(kappa == 0, theta(t, 0), theta(t, 1))
    W = rot * np.exp(1j * th) * L
    Z = W.real.copy()
    imrat = max(imrat, float(np.max(np.abs(W.imag[1:]) / (np.abs(W[1:]) + 1e-3))))
    if i == pad: z0row = Z.copy()
    buf.append(Z); rowid.append(i)
    if len(buf) == 8:                         # rows r-3 .. r+4 with r = local index 3
        za, zb = buf[3], buf[4]
        cols = np.nonzero(za * zb < 0)[0]; cols = cols[cols >= 1]
        if len(cols):
            vals = np.stack([buf[d][cols] for d in range(8)], axis=1)
            u = newton(vals @ V7.T, za[cols] / (za[cols] - zb[cols]))
            u5 = newton(vals[:, 1:7] @ V5.T, u.copy())
            gam = ts[rowid[3]] + u * h
            KZ.append(cols.astype(np.int32)); GZ.append(gam); ERR.append((np.abs(u5 - u) * h).astype(np.float32))
    if i % 500 == 0:
        print(f'row {i}/{len(ts)} t = {t:.3f}  {time.time() - t0:.0f} s', flush=True)
kz = np.concatenate(KZ); gz = np.concatenate(GZ); er = np.concatenate(ERR)
print(f'grid done: {len(ts)} rows, {time.time() - t0:.0f} s; max |Im W|/(|W| + 1e-3) = {imrat:.2e}; '
      f'interpolation difference deg 7 vs 5: max {er.max():.2e}, median {np.median(er):.2e}', flush=True)
ok = (gz > 0) & (gz <= T)
kz, gz, er = kz[ok], gz[ok], er[ok]
found = np.bincount(kz, minlength=q - 1)
kbar = (-k) % (q - 1)
zero_at_0 = np.abs(z0row) < 1e-9
# argument principle at t = T, streamed over sigma
sig = np.concatenate([np.linspace(2.0, 0.6, 141), np.linspace(0.6, 0.5, 201)[1:]])
prev = np.angle(L_all(sig[0] + 1j * T)); A = prev.copy(); maxstep = 0.0
for s_ in sig[1:]:
    cur = np.angle(L_all(s_ + 1j * T))
    d = (cur - prev + math.pi) % (2 * math.pi) - math.pi
    maxstep = max(maxstep, float(np.abs(d[1:]).max())); A += d; prev = cur
NAP = (2 * np.where(kappa == 0, theta(T, 0), theta(T, 1)) + A + A[kbar]) / math.pi
NAPr = np.rint(NAP); dev = np.abs(NAP - NAPr)[1:].max()
count = found + found[kbar] + zero_at_0
mism = np.nonzero((NAPr != count) & (k >= 1))[0]
print(f'characters with |Z_k(0)| < 1e-9: {int(zero_at_0[1:].sum())}; argument principle: largest step {maxstep:.3f} rad, '
      f'largest distance from an integer {dev:.2e}; characters whose count differs: {len(mism)} {mism[:10].tolist()}', flush=True)
print(f'zeros with 0 < gamma <= T: {len(gz)}; with |gamma| <= T over all chi: {int(count[1:].sum())}; '
      f'N(T,chi) from {int(count[1:].min())} to {int(count[1:].max())}; {time.time() - t0:.0f} s', flush=True)
np.savez(out, q=q, g=g, T=T, h=h, k=kz, gamma=gz, interr=er, count=count, NAP=NAP, pw=pw, mism=mism)
