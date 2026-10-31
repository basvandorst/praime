# Puzzle 374: gaps between self numbers. S(10q + r) = B_q + 2r, B_q = 10q + sod(q); B_(q+1) - B_q = 11 - 9 t(q).
import numpy as np
from collections import Counter
N = 10**7
sod = np.zeros(N + 200, np.int64)
for i in range(1, N + 200): sod[i] = sod[i // 10] + i % 10
gen = np.zeros(N + 400, bool)
j = np.arange(N + 200); S = j + sod[j]; gen[S[S < N + 400]] = True
selfs = np.nonzero(~gen[1:N])[0] + 1
gaps = np.diff(selfs)
c = Counter(gaps.tolist())
tot = len(gaps)
print("self numbers below 1e7:", len(selfs), " gaps:", tot)
print("most common gaps:", [(g, n, f"{n / tot:.4f}") for g, n in c.most_common(8)])
print("gap 11 share:", f"{c[11] / tot:.4f}", " all gaps 11 or 2 mod 13 below 1e7:", all(g == 11 or g % 13 == 2 for g in c))
# relate each gap to the block structure: self number n = B_q + 9 in a 'normal' stretch
B = lambda q: 10 * q + sod[q]
def t(q):
    k = 0
    while q % 10 == 9: q //= 10; k += 1
    return k
# where the gap is not 11: the q of the block below the self number, and t of the previous block
odd = Counter()
for a, g in zip(selfs[:-1], gaps):
    if g == 11: continue
    q = (a - 9) // 10                 # rough block index
    # find the block q with B_q <= a < B_q + 20 and B_q + 9 closest
    tq = max(t(qq) for qq in range(max(q - 3, 0), q + 2))
    odd[(int(g), tq)] += 1
print("non-11 gaps by (gap, max trailing nines nearby):", sorted(odd.items())[:20])
