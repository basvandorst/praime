"""All central values L(1/2, chi) mod 395,137 from one FFT; their extremes, their link with the lowest zero and with
the small primes (Selberg: log|L(1/2,chi)| is about sum_p Re chi(p)/sqrt p)."""
import math, numpy as np
from scipy.special import bernoulli
D = np.load('lfam395137f.npz'); q = int(D['q']); pw = D['pw']; kz = D['k'].astype(np.int64); gz = D['gamma']
N = q - 1; ind = np.empty(q, dtype=np.int64); ind[pw] = np.arange(N)
NT_, M_ = 30, 9; alpha = pw / q
LOG = np.log(np.arange(NT_)[None, :] + alpha[:, None]); Na = NT_ + alpha; LNa = np.log(Na)
B = bernoulli(2 * M_); coef = [float(B[2 * j]) / math.factorial(2 * j) for j in range(1, M_ + 1)]
s = 0.5 + 0j
main = np.exp(-s * LOG).sum(axis=1); Nas = np.exp(-s * LNa)
tail = Nas * Na / (s - 1) + 0.5 * Nas; rising = s; powr = Nas / Na
for j in range(1, M_ + 1):
    tail = tail + coef[j - 1] * rising * powr; rising = rising * (s + 2 * j - 1) * (s + 2 * j); powr = powr / (Na * Na)
Lc = np.exp(-s * math.log(q)) * N * np.fft.ifft(main + tail)       # L(1/2, chi_k), k = 0 .. N-1
ks = np.arange(1, N); A = np.abs(Lc[ks])
low = np.full(N, np.inf); np.minimum.at(low, kz, gz); lowest = np.minimum(low[ks], low[(-ks) % N])
primes = [p for p in range(2, 1000) if all(p % r for r in range(2, int(p ** 0.5) + 1))]
sel = sum(np.cos(2 * math.pi * ((ks * ind[p]) % N) / N) / math.sqrt(p) for p in primes)
la = np.log(A)
out = [f'central values of the {len(ks):,} non-principal characters: |L(1/2)| min {A.min():.3e} (chi_{int(ks[np.argmin(A)])}), max {A.max():.2f} (chi_{int(ks[np.argmax(A)])}), median {np.median(A):.3f}',
       f'log|L(1/2)|: mean {la.mean():+.3f}, variance {la.var():.3f} (Selberg: variance about (1/2) ln ln q = {0.5 * math.log(math.log(q)):.3f})',
       f'correlation of log|L(1/2)| with sum_(p<1000) Re chi(p)/sqrt p: {np.corrcoef(la, sel)[0, 1]:+.3f}; with the lowest zero: {np.corrcoef(la, lowest)[0, 1]:+.3f}; with log of the lowest zero: {np.corrcoef(la, np.log(lowest))[0, 1]:+.3f}',
       f'the 5 smallest |L(1/2)|: ' + ', '.join(f'chi_{int(ks[i])}: {A[i]:.2e} (lowest zero {lowest[i]:.2e})' for i in np.argsort(A)[:5]),
       f'the 5 largest |L(1/2)|: ' + ', '.join(f'chi_{int(ks[i])}: {A[i]:.2f} (lowest zero {lowest[i]:.3f})' for i in np.argsort(A)[::-1][:5]),
       f'moments: mean |L|^2 = {np.mean(A ** 2):.4f}, mean |L|^4 = {np.mean(A ** 4):.3f} (ln q = {math.log(q):.3f})']
# law: lowest zero against |L(1/2)| for small central values (a simple zero near the centre gives |L(1/2)| ~ |L'| gamma_1)
s_ = A < 1e-2
out.append(f'characters with |L(1/2)| < 0.01: {int(s_.sum())}; for them corr(log|L(1/2)|, log lowest zero) = {np.corrcoef(la[s_], np.log(lowest[s_]))[0, 1]:+.3f}')
print('\n'.join(out)); open('probeA3.txt', 'w').write('\n'.join(out) + '\n')
