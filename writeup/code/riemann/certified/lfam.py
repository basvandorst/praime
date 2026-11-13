"""Zeros of all Dirichlet L-functions modulo a prime q on 0 <= t <= T, with a completeness check.

L(s, chi_k) = q^-s sum_{a=1}^{q-1} chi_k(a) zeta(s, a/q), chi_k(g^j) = exp(2 pi i j k/(q-1)), g a primitive root.
The Hurwitz zeta values come from Euler-Maclaurin (N direct terms, M Bernoulli terms) for all a at once;
the sum over a is one FFT in the discrete-log order. Z_k(t) = Re(eps_k^(-1/2) e^{i theta_k(t)} L(1/2+it, chi_k)) is
real (eps_k = tau(chi_k)/(i^kappa sqrt q)); its sign changes on a grid of step h give the zeros, refined by
degree-7 interpolation. Completeness: for each chi the argument principle gives
  N(T, chi) = #{|gamma| <= T} = (2 theta(T) + arg L(1/2+iT, chi) + arg L(1/2+iT, chibar)) / pi,
with arg L tracked from sigma = 2 down to 1/2; it must equal the number of sign changes found for chi on (0, T]
plus those for chibar. Usage: lfam.py q T h out.npz
"""
import sys, math, time
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
kk = np.arange(NT_)
LOG = np.log(kk[None, :] + alpha[:, None])           # (q-1, N)
Na = NT_ + alpha; LNa = np.log(Na)
B = bernoulli(2 * M_)                                 # B[2j]
coef = [float(B[2 * j]) / math.factorial(2 * j) for j in range(1, M_ + 1)]

def hurwitz_all(s):
    """zeta(s, alpha_j) for all j (complex s, |s| up to ~20)."""
    main = np.exp(-s * LOG).sum(axis=1)
    Nas = np.exp(-s * LNa)                            # (N+alpha)^-s
    tail = Nas * Na / (s - 1) + 0.5 * Nas
    rising = s; powr = Nas / Na                       # (s)_{2j-1}, (N+alpha)^(-s-2j+1)
    for j in range(1, M_ + 1):
        tail = tail + coef[j - 1] * rising * powr
        rising = rising * (s + 2 * j - 1) * (s + 2 * j)
        powr = powr / (Na * Na)
    return main + tail

def L_all(s):
    H = hurwitz_all(s)
    return np.exp(-s * math.log(q)) * (q - 1) * np.fft.ifft(H)

k = np.arange(q - 1)
kappa = k % 2
tau = (q - 1) * np.fft.ifft(np.exp(2j * math.pi * pw / q))
eps = tau / ((1j) ** kappa * math.sqrt(q))
print('|eps| range (k >= 1):', np.abs(eps[1:]).min(), np.abs(eps[1:]).max(), flush=True)
rot = np.exp(-0.5j * np.angle(eps))

def theta(t, kap):
    return 0.5 * t * math.log(q / math.pi) + np.imag(loggamma((0.5 + kap + 1j * t) / 2))

t0 = time.time()
pad = 8
ts = (np.arange(int(math.ceil(T / h)) + 2 * pad + 1) - pad) * h
Z = np.empty((len(ts), q - 1))
imrat = 0.0
for i, t in enumerate(ts):
    L = L_all(0.5 + 1j * t)
    th = np.where(kappa == 0, theta(t, 0), theta(t, 1))
    W = rot * np.exp(1j * th) * L
    Z[i] = W.real
    imrat = max(imrat, float(np.max(np.abs(W.imag[1:]) / (np.abs(W[1:]) + 1e-3))))
    if i % 500 == 0:
        print(f'row {i}/{len(ts)} t = {t:.3f}  {time.time() - t0:.0f} s', flush=True)
print(f'grid done: {len(ts)} rows, {time.time() - t0:.0f} s; max |Im W|/(|W| + 1e-3) = {imrat:.2e}', flush=True)

# zeros: sign changes, refined with an 8-point interpolant (rows i-3 .. i+4, local u in -3 .. 4)
U = np.arange(-3, 5, dtype=float)
Vinv = np.linalg.inv(np.vander(U, 8, increasing=True))
rows, cols = np.nonzero(Z[:-1] * Z[1:] < 0)
sel = (rows >= 3) & (rows + 4 < len(ts)) & (cols >= 1)
rows, cols = rows[sel], cols[sel]
vals = np.stack([Z[rows + d, cols] for d in range(-3, 5)], axis=1)
C = vals @ Vinv.T                                     # coefficients in u (increasing powers)
z0, z1 = Z[rows, cols], Z[rows + 1, cols]
u = z0 / (z0 - z1)                                    # linear guess in [0, 1]
for _ in range(30):
    p = np.zeros_like(u); dp = np.zeros_like(u)
    for c_ in C.T[::-1]:
        dp = dp * u + p; p = p * u + c_
    u = np.clip(u - p / dp, 0, 1)
gam = ts[rows] + u * h
# error estimate: degree-5 interpolant on rows i-2 .. i+3
U5 = np.arange(-2, 4, dtype=float); V5 = np.linalg.inv(np.vander(U5, 6, increasing=True))
C5 = np.stack([Z[rows + d, cols] for d in range(-2, 4)], axis=1) @ V5.T
u5 = u.copy()
for _ in range(30):
    p = np.zeros_like(u5); dp = np.zeros_like(u5)
    for c_ in C5.T[::-1]:
        dp = dp * u5 + p; p = p * u5 + c_
    u5 = np.clip(u5 - p / dp, 0, 1)
interr = np.abs(u5 - u) * h
print(f'{len(gam)} sign changes; interpolation difference deg 7 vs 5: max {interr.max():.2e}, median {np.median(interr):.2e}', flush=True)
ok = (gam > 0) & (gam <= T)
gk, gg = cols[ok], gam[ok]
found = np.bincount(gk, minlength=q - 1)
kbar = (-k) % (q - 1)
zero_at_0 = np.abs(Z[pad]) < 1e-9
print('characters with |Z_k(0)| < 1e-9:', int(np.sum(zero_at_0[1:])))

# argument principle at t = T for all k
sig = np.concatenate([np.linspace(2.0, 0.6, 141), np.linspace(0.6, 0.5, 201)[1:]])
args = []
for s_ in sig:
    args.append(np.angle(L_all(s_ + 1j * T)))
args = np.array(args)
d = np.diff(args, axis=0); d = (d + math.pi) % (2 * math.pi) - math.pi
maxstep = np.abs(d[:, 1:]).max()
A = args[0] + d.sum(axis=0)
NAP = (2 * np.where(kappa == 0, theta(T, 0), theta(T, 1)) + A + A[kbar]) / math.pi
NAPr = np.rint(NAP)
dev = np.abs(NAP - NAPr)[1:].max()
count = found + found[kbar] + zero_at_0
mism = np.nonzero((NAPr != count) & (k >= 1))[0]
print(f'argument principle: largest step {maxstep:.3f} rad, largest distance from an integer {dev:.2e}; '
      f'characters whose count differs: {len(mism)} {mism[:10]}', flush=True)
print(f'zeros with 0 < gamma <= T: {len(gg)}; with |gamma| <= T over all chi: {int(count[1:].sum())}; '
      f'N(T,chi) found from {int(count[1:].min())} to {int(count[1:].max())}; {time.time() - t0:.0f} s', flush=True)
np.savez(out, q=q, g=g, T=T, h=h, k=gk, gamma=gg, interr=interr[ok], count=count, NAP=NAP, pw=pw)
