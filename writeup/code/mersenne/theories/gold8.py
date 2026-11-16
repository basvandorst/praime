"""Theory 9, second modulus (8): every odd square is 1 mod 8, so the class 1 mod 8 carries the prime squares s^2 (s >= 3).
  n = 0 mod 8: G(1,7)+G(7,1) - G(3,5)-G(5,3) ~ -S(n);  n = 2: G(1,1) - G(5,5) ~ -S(n);  n = 4: G(1,3)+G(3,1) - G(5,7)-G(7,5) ~ -S(n);
  n = 6: G(1,5)+G(5,1) - G(3,3)-G(7,7) ~ -S(n);  S(n) = #{prime s >= 3 : n - s^2 prime}.  (Otherwise as gold.py.)
Theory 9: Chebyshev's bias in Goldbach's comet. All prime squares s^2 (s >= 5) are 1 mod 12. If the prime-power count
(primes plus half a point at each prime square) is unbiased in additive problems too, then for even n the ordered
representations n = p + q that use the class 1 mod 12 lose half a pair for each prime s with n - s^2 prime, on each side:
  n = 0 mod 12: G(1,11)+G(11,1) - G(5,7)-G(7,5)  ~ -S(n)
  n = 2 mod 12: G(1,1) - G(7,7)                  ~ -S(n)
  n = 6 mod 12: G(1,5)+G(5,1) - G(7,11)-G(11,7)  ~ -S(n)
  n = 10 mod 12: G(5,5) - G(11,11)               ~ 0 (control: no class 1)
S(n) = #{prime s >= 5 : s^2 < n, n - s^2 prime}. usage: gold.py N"""
import sys, numpy as np
N = int(float(sys.argv[1]))
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
n = np.arange(N + 1)
L = 1 << int(np.ceil(np.log2(2 * N + 2)))
F = {a: np.fft.rfft((s & (n % 8 == a)).astype(np.float64), L) for a in (1, 3, 5, 7)}
def G(a, b): return np.rint(np.fft.irfft(F[a] * F[b], L)[:N + 1]).astype(np.int64)
S = np.zeros(N + 1, np.int64)
for q in np.nonzero(s)[0]:
    if q < 3: continue
    if q * q >= N: break
    S[q * q + 1:] += s[1:N + 1 - q * q]
lo = N // 4
D = {0: G(1, 7) + G(7, 1) - G(3, 5) - G(5, 3), 2: G(1, 1) - G(5, 5), 4: G(1, 3) + G(3, 1) - G(5, 7) - G(7, 5), 6: G(1, 5) + G(5, 1) - G(3, 3) - G(7, 7)}
tot = {0: G(3, 5) + G(5, 3), 2: G(5, 5), 4: G(5, 7) + G(7, 5), 6: G(3, 3) + G(7, 7)}
chk = np.arange(lo, N + 1)
chk = chk[(chk % 2 == 0)][:2000]
brute_ok = True
P = np.nonzero(s)[0]
for m in chk[:200]:
    c = m % 8
    if c not in D: continue
    pp = P[P < m]; qq = m - pp; ok = s[qq]
    cls = lambda a, b: int(np.sum(ok & (pp % 8 == a) & (qq % 8 == b)))
    want = {0: cls(1, 7) + cls(7, 1) - cls(3, 5) - cls(5, 3), 2: cls(1, 1) - cls(5, 5), 4: cls(1, 3) + cls(3, 1) - cls(5, 7) - cls(7, 5), 6: cls(1, 5) + cls(5, 1) - cls(3, 3) - cls(7, 7)}[c]
    brute_ok &= want == D[c][m]
print(f'N = {N:,}; n in [{lo:,}, {N:,}]; FFT counts equal a direct count on 200 values of n: {brute_ok}')
print('Identity (proved): Delta(n) = Delta_Pi(n) - S(n), where Delta_Pi counts the class 1 with half a point at each prime square. Theory: Delta_Pi has mean 0 (no bias left).')
print('range                n mod 8    count   mean Delta   mean -S(n)   mean Delta_Pi = Delta + S   sd of that mean   corr(Delta_Pi/G, S/G)')
for a, b in ((N // 16, N // 8), (N // 8, N // 4), (N // 4, N // 2), (N // 2, N + 1)):
    for c in (0, 2, 4, 6):
        m = (n >= a) & (n < b) & (n % 8 == c)
        d = D[c][m].astype(float); x = S[m].astype(float); g = tot[c][m].astype(float)
        dp = d + x
        cc = np.corrcoef(dp / g, x / g)[0, 1]
        print(f'[{a:>9,}, {b:>9,}) {c:>6} {m.sum():>9,} {d.mean():>11.2f} {-x.mean():>12.2f} {dp.mean():>23.2f} {dp.std() / np.sqrt(len(dp)):>15.2f} {cc:>16.3f}')
