"""Theory 9b: the rest of the Goldbach bias is the prime race. To first order in the class fluctuations,
  n = 0 mod 12: Delta ~ 2 w [(pi_1 - pi_5) + (pi_11 - pi_7)](n)     n = 2: Delta ~ 2 w (pi_1 - pi_7)(n)
  n = 6:        Delta ~ 2 w [(pi_1 - pi_7) + (pi_5 - pi_11)](n)     n = 10: Delta ~ 2 w (pi_5 - pi_11)(n)
(pi_a(n) = primes <= n that are a mod 12, w = the chance that the partner is prime). The prime squares are the part of the race
that is known exactly: each prime square s^2 <= n (s >= 5) takes half a prime from class 1, i.e. -Q(n)/2 in pi_1 - pi_7, and in
Goldbach it gives exactly -S(n). So w = S(n)/Q(n), and with no parameter
  Delta(n) ~ -S(n) * D(n) / (Q(n)/2),  D = the lead of the other classes over class 1 as above (n = 10: (pi_11 - pi_5), w from n - 2).
Q(n) = #{prime s >= 5: s^2 <= n}. Compared in 40 bins of log n over [N/16, N]. usage: goldrace.py N [M]
General form, M = 12 or 8 (every prime square from s >= 5, resp. s >= 3, is 1 mod M): for pair types (a, b) with signs,
Delta ~ (S/Q) * sum of sign * (pi_a + pi_b), where pi_1 appears twice with sign +."""
import sys, numpy as np
N = int(float(sys.argv[1])); M = int(sys.argv[2]) if len(sys.argv) > 2 else 12
CL = (1, 5, 7, 11) if M == 12 else (1, 3, 5, 7); S0 = 5 if M == 12 else 3
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
n = np.arange(N + 1)
L = 1 << int(np.ceil(np.log2(2 * N + 2)))
F = {a: np.fft.rfft((s & (n % M == a)).astype(np.float64), L) for a in CL}
def G(a, b): return np.rint(np.fft.irfft(F[a] * F[b], L)[:N + 1]).astype(np.int64)
pi = {a: np.cumsum(s & (n % M == a)) for a in CL}
S = np.zeros(N + 1, np.int64); Q = np.zeros(N + 1, np.int64)
for q in np.nonzero(s)[0]:
    if q < S0: continue
    if q * q > N: break
    S[q * q + 1:] += s[1:N + 1 - q * q]
    Q[q * q:] += 1
TYPES = {12: {0: [((1, 11), 1), ((11, 1), 1), ((5, 7), -1), ((7, 5), -1)], 2: [((1, 1), 1), ((7, 7), -1)],
              6: [((1, 5), 1), ((5, 1), 1), ((7, 11), -1), ((11, 7), -1)], 10: [((5, 5), 1), ((11, 11), -1)]},
         8: {0: [((1, 7), 1), ((7, 1), 1), ((3, 5), -1), ((5, 3), -1)], 2: [((1, 1), 1), ((5, 5), -1)],
             4: [((1, 3), 1), ((3, 1), 1), ((5, 7), -1), ((7, 5), -1)], 6: [((1, 5), 1), ((5, 1), 1), ((3, 3), -1), ((7, 7), -1)]}}[M]
D = {c: sum(sg * G(a, b) for (a, b), sg in T) for c, T in TYPES.items()}
SIG = {c: sum(sg * (pi[a] + pi[b]) for (a, b), sg in T) for c, T in TYPES.items()}
ctrl = [c for c, T in TYPES.items() if all(1 not in ab for ab, _ in T)]
Sw = S.astype(float).copy(); Sw[2:] = S[:-2]
with np.errstate(divide='ignore', invalid='ignore'):
    pred = {c: (Sw if c in ctrl else S) * SIG[c] / Q for c in TYPES}
edges = np.exp(np.linspace(np.log(N / 16), np.log(N), 41))
print(f'N = {N:,}; bins of log n from {N // 16:,}; prediction Delta ~ -S(n) * D(n) / (Q(n)/2), no parameter')
print(f' n mod {M}  corr over bins (race)  corr (squares only, -S)  mean |error| race  mean |error| squares only  sd of a bin mean')
rows = []
for c in TYPES:
    bd, bp, bs, se = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (n >= a) & (n < b) & (n % M == c) & (Q > 0)
        d = D[c][m].astype(float)
        bd.append(d.mean()); bp.append(pred[c][m].mean()); bs.append(-S[m].mean()); se.append(d.std() / np.sqrt(len(d)))
    bd, bp, bs, se = map(np.array, (bd, bp, bs, se))
    cs = np.corrcoef(bd, bs)[0, 1] if c not in ctrl else float('nan')
    print(f'{c:>9} {np.corrcoef(bd, bp)[0, 1]:>22.3f} {cs:>24.3f} {np.mean(np.abs(bd - bp)):>18.2f} {np.mean(np.abs(bd - bs)):>26.2f} {np.mean(se):>17.2f}')
    rows.append((c, bd, bp, bs))
print(f'\nper bin (n = 0 mod {M}): lower edge, mean Delta, race prediction, squares only')
c, bd, bp, bs = rows[0]
for a, x, y, z in zip(edges[:-1], bd, bp, bs):
    print(f'  {a:>12,.0f} {x:>9.2f} {y:>9.2f} {z:>9.2f}')
