"""Puzzle 9: the numbers behind the heuristic model in writeup/9.md.

Sieves the primes below 10^8 (k <= 5761455) and measures
- the asymptotic S(k) ~ (k^2/2)(ln k + ln ln k - 3/2 + (ln ln k - 5/2)/ln k),
- the local factor E[rho_q(S(k) mod q)] for q = 8, 3, 5, 7, 840 (1 = no bias),
- the expected number of squares sum 1/(2 sqrt S(k)) against the closed form sqrt(2 ln K),
- the expected number of m-th powers for m = 3, 5, 7 and the tail beyond K,
then prints the forecast for squares beyond Resta's k = 126789311423.
Run from the repo root: .venv/bin/python writeup/code/p9_model.py   (seconds, ~1 GB)
"""
import math
import numpy as np

N = 10**8
s = np.ones(N // 2, dtype=bool)
s[0] = False
for i in range(3, int(N**0.5) + 1, 2):
    if s[i // 2]:
        s[i * i // 2::i] = False
p = np.concatenate(([2], 2 * np.nonzero(s)[0] + 1)).astype(np.int64)
S = np.cumsum(p)
T = S[1:] - 2
K = len(p)
print("K =", K, " S(K) =", int(S[-1]))

print("\nasymptotic: S(k) / ((k^2/2)(ln k + lnln k - 3/2 + (lnln k - 5/2)/ln k))")
for k in (2474, 709838, K):
    L = math.log(k)
    LL = math.log(L)
    print("  k = %-8d ratio %.6f" % (k, S[k - 1] / (k * k / 2 * (L + LL - 1.5 + (LL - 2.5) / L))))


def rho(q):
    r = np.zeros(q)
    for y in range(q):
        r[y * y % q] += 1
    return r


print("\nlocal factor E[rho_q] (number of square roots of S mod q, averaged over k)")
for name, A in (("S", S), ("T", T)):
    print("  %s: " % name + ", ".join("q=%d %.4f" % (q, rho(q)[(A % q).astype(np.int64)].mean())
                                     for q in (8, 3, 5, 7, 840)))

Sf, Tf = S.astype(np.float64), T.astype(np.float64)
print("\nexpected squares sum_{k<=K} 1/(2 sqrt S(k)) vs closed form sqrt(2 ln K) - sqrt(2 ln 2)")
for hi in (10, 10**4, 10**6, K):
    print("  K = %-8d sum %.3f   closed form %.3f" % (hi, (0.5 / np.sqrt(Sf[:hi])).sum(),
                                                     math.sqrt(2 * math.log(hi)) - math.sqrt(2 * math.log(2))))
print("  T: sum %.3f" % (0.5 / np.sqrt(Tf)).sum())

print("\nexpected m-th powers: sum_{k<=K} (1/m) S^(1/m - 1), and the tail beyond K")
for m in (3, 5, 7):
    L = math.log(K)
    tail = 2 ** ((m - 1) / m) * K ** (-(m - 2) / m) / ((m - 2) * L ** ((m - 1) / m))
    print("  m = %d: S %.4f  T %.4f  tail beyond K %.2e" % (m, (Sf ** (1 / m - 1)).sum() / m,
                                                          (Tf ** (1 / m - 1)).sum() / m, tail))
for B in (1.9e8, 3.6e12):
    L = math.log(B)
    print("  cube tail beyond k = %.1e: %.2e" % (B, 2 ** (2 / 3) * B ** (-1 / 3) / L ** (2 / 3)))

print("\nforecast beyond Resta (k1 = 126789311423): E = sqrt(2 ln K) - sqrt(2 ln k1), P(>=1) = 1 - e^-E")
k1 = 126789311423
for Kf in (1e12, 1e13, 2.3e13, 1e14, 1e15, 1e16, 1e18):
    E = math.sqrt(2 * math.log(Kf)) - math.sqrt(2 * math.log(k1))
    pK = Kf * (math.log(Kf) + math.log(math.log(Kf)) - 1)
    print("  K = %.1e  E = %.3f  P = %.3f  p_K ~ %.2e" % (Kf, E, 1 - math.exp(-E), pK))
