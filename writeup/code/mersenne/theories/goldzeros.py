"""Theory 9c: Goldbach's class strands from the zeros alone. In characters mod 12 the race of theory 9b is
  n = 0: SIG = 2 sum chi_12(p);  n = 2: sum (chi_-4 + chi_12)(p);  n = 6: 2 sum chi_-4(p);  n = 10: sum (chi_-4 - chi_12)(p)
(p >= 5). For a real character, sum_{p <= n} chi(p) = Pi(n, chi) - (1/2) #{p <= sqrt n, chi(p^2) = 1} - ..., and the explicit
formula gives Pi(n, chi) ~ -Z_chi(n), Z_chi(n) = sum over the zeros rho of L(s, chi) of Li(n^rho) (here 289 zeros of L(s, chi_-4)
and 359 of L(s, chi_12) below height 400, from PARI's lfunzeros; Li(x^rho) ~ x^rho/(rho ln x) (1 + 1/(rho ln x) + 2/(rho ln x)^2)).
So, with the squares exact (Q = #{prime s >= 5, s^2 <= n}; chi_-4 also has 3^2):
  n = 0: Delta ~ -S (1 + 2 Z_12/Q)            n = 2: Delta ~ -S (1 + (Z_-4 + Z_12 - 1/2)/Q)
  n = 6: Delta ~ -S (1 + (2 Z_-4 - 1)/Q)       n = 10: Delta ~ S(n-2) (Z_12 - Z_-4 + 1/2)/Q
No prime count enters the race part; the zeros replace it. usage: goldzeros.py N"""
import sys, numpy as np
N = int(float(sys.argv[1]))
Z4 = np.loadtxt('zeros_m4.txt'); Z12 = np.loadtxt('zeros_12.txt')
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
n = np.arange(N + 1)
L = 1 << int(np.ceil(np.log2(2 * N + 2)))
F = {a: np.fft.rfft((s & (n % 12 == a)).astype(np.float64), L) for a in (1, 5, 7, 11)}
def G(a, b): return np.rint(np.fft.irfft(F[a] * F[b], L)[:N + 1]).astype(np.int64)
S = np.zeros(N + 1, np.int64); Q = np.zeros(N + 1, np.int64)
for q in np.nonzero(s)[0]:
    if q < 5: continue
    if q * q > N: break
    S[q * q + 1:] += s[1:N + 1 - q * q]
    Q[q * q:] += 1
D = {0: G(1, 11) + G(11, 1) - G(5, 7) - G(7, 5), 2: G(1, 1) - G(7, 7), 6: G(1, 5) + G(5, 1) - G(7, 11) - G(11, 7), 10: G(5, 5) - G(11, 11)}
pi = {a: np.cumsum(s & (n % 12 == a)) for a in (1, 5, 7, 11)}
def Zsum(x, g):
    x = x.astype(np.float64); lx = np.log(x)
    out = np.zeros(len(x))
    rho = 0.5 + 1j * g
    for i in range(0, len(x), 20000):
        l = lx[i:i + 20000, None]
        u = rho[None, :] * l
        t = np.exp(u) / u * (1 + 1 / u + 2 / u ** 2)
        out[i:i + 20000] = 2 * t.real.sum(1)
    return out
edges = np.exp(np.linspace(np.log(N / 16), np.log(N), 41))
print(f'N = {N:,}; zeros below 400: {len(Z4)} of L(s, chi_-4) (first {Z4[0]:.4f}), {len(Z12)} of L(s, chi_12) (first {Z12[0]:.4f})')
print(' n mod 12   corr: zeros   corr: measured race   corr: squares only   mean |err|: zeros   race   squares only')
out = {}
for c in (0, 2, 6, 10):
    m = (n >= edges[0]) & (n % 12 == c) & (Q > 0)
    x = n[m]; z4 = Zsum(x, Z4); z12 = Zsum(x, Z12); q = Q[m].astype(float); sv = S[m].astype(float)
    if c == 0: pz = -sv * (1 + 2 * z12 / q); sig = 2 * ((pi[1] + pi[11] - pi[5] - pi[7])[m])
    if c == 2: pz = -sv * (1 + (z4 + z12 - 0.5) / q); sig = 2 * (pi[1] - pi[7])[m]
    if c == 6: pz = -sv * (1 + (2 * z4 - 1) / q); sig = 2 * ((pi[1] - pi[7]) + (pi[5] - pi[11]))[m]
    if c == 10:
        sv = S[x - 2].astype(float); pz = sv * (z12 - z4 + 0.5) / q; sig = 2 * (pi[5] - pi[11])[m]
    pr = sv * sig / q
    d = D[c][m].astype(float)
    bd, bz, br, bs = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (x >= a) & (x < b)
        bd.append(d[k].mean()); bz.append(pz[k].mean()); br.append(pr[k].mean()); bs.append(-S[x[k]].mean() if c != 10 else 0.0)
    bd, bz, br, bs = map(np.array, (bd, bz, br, bs))
    cs = np.corrcoef(bd, bs)[0, 1] if c != 10 else float('nan')
    print(f'{c:>9} {np.corrcoef(bd, bz)[0, 1]:>13.3f} {np.corrcoef(bd, br)[0, 1]:>21.3f} {cs:>20.3f} {np.mean(np.abs(bd - bz)):>19.2f} {np.mean(np.abs(bd - br)):>6.2f} {np.mean(np.abs(bd - bs)):>14.2f}', flush=True)
    out[c] = (bd, bz, br, bs)
print('\nper bin (n = 0 mod 12): lower edge, mean Delta, from the zeros, from the measured race, squares only')
bd, bz, br, bs = out[0]
for a, w, xx, y, z in zip(edges[:-1], bd, bz, br, bs):
    print(f'  {a:>12,.0f} {w:>9.2f} {xx:>9.2f} {y:>9.2f} {z:>9.2f}')
