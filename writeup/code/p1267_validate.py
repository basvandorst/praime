"""Puzzle 1267: how many emirps are there?  Counts the emirps (prime whose digit reversal is a different prime) with
k digits for k <= 8 by a sieve to 10^8, and compares with the model
    E(k) = sum over k-digit primes p with leading digit in {1, 3, 7, 9}, p not a palindrome, of (33/8) / ln(rev p),
where 33/8 = 2 * (3/2) * (5/4) * (11/10) accounts for rev p being odd, prime to 3 (same digit sum as p), prime to 5
(its last digit is the first digit of p) and prime to 11 (alternating digit sum is +-that of p).
Then the closed form of the model with the logarithmic integral is compared with the page's counts for k <= 12.
"""
import numpy as np, math
from mpmath import li, mpf, log

N = 10**8
s = np.ones(N + 1, dtype=bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]:
        s[i*i::i] = False
cnt, model, pal = {}, {}, {}
for p in np.flatnonzero(s).tolist():
    st = str(p); r = int(st[::-1]); k = len(st)
    if r == p:
        pal[k] = pal.get(k, 0) + 1
        continue
    if st[0] in '1379':
        model[k] = model.get(k, 0.0) + 4.125 / math.log(r)
    if s[r]:
        cnt[k] = cnt.get(k, 0) + 1
page = {2: 8, 3: 28, 4: 204, 5: 1406, 6: 9538, 7: 70474, 8: 535578, 9: 4192024, 10: 33619380, 11: 274890230, 12: 2294771254}
assert all(cnt[k] == page[k] for k in range(2, 9)), cnt
print("emirps with k digits (sieve to 10^8):", {k: cnt[k] for k in sorted(cnt)})
print("palindromic primes with k digits:", pal)
print("exact model sum:", {k: round(model[k]) for k in sorted(model)})
print("actual / exact model:", {k: round(cnt[k] / model[k], 4) for k in sorted(cnt)})

# closed form: (33/8) * sum_{d in 1,3,7,9} [li((d+1) 10^(k-1)) - li(d 10^(k-1))] * mean_{d'} 1/ln(d' 10^(k-1) + 5 10^(k-2))
print("closed form vs page counts:")
for k in range(2, 13):
    B = mpf(10)**(k - 1)
    lead = sum(li((d + 1) * B) - li(d * B) for d in (1, 3, 7, 9))
    inv = sum(1 / log(dp * B + B / 2) for dp in (1, 3, 7, 9)) / 4
    E = mpf(33) / 8 * lead * inv
    print(f"  k={k:2d}: model {float(E):14.0f}  page {page[k]:12d}  ratio {page[k] / float(E):.4f}")
# asymptotic constant: emirps(x) ~ (33/8)(4/9) x / ln^2 x = (11/6) x / ln^2 x
x = mpf(10)**12; print("cumulative to 10^12:", sum(page.values()), " (11/6) x/ln^2 x =", float(mpf(11)/6 * x / log(x)**2))
