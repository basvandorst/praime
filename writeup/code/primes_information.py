import sys, time
import numpy as np
from math import log, log2, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve

X = 10**8
t0 = time.time()
isp = sieve(X)
pr = np.flatnonzero(isp)
g = np.diff(pr)


def H(vals):
    _, c = np.unique(vals, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def Hcond(keys, vals):
    kv = keys.astype(np.int64) * (vals.max() + 1) + vals
    return H(kv) - H(keys)


print("E-F1 information in the gap sequence per scale (gaps of the primes in [x, 2x]), bits per gap:")
print("     x      H(gap)   H(gap | p mod 30)   H(gap | p mod 210)   H(gap | last gap)   log2(ln x)+const fit")
rows = []
for x in (10**4, 10**5, 10**6, 10**7, 5 * 10**7):
    a, b = np.searchsorted(pr, [x, 2 * x])
    gg = g[a:b - 1]; pp = pr[a:b - 1]
    cap = 6 * 200
    gc = np.minimum(gg, cap) // 2
    h0 = H(gc); h30 = Hcond(pp % 30, gc); h210 = Hcond(pp % 210, gc); hl = Hcond(np.concatenate(([0], gc[:-1])), gc)
    rows.append((log2(log(x)), h0))
    print(f"  {x:9d}   {h0:.3f}      {h30:.3f}              {h210:.3f}               {hl:.3f}")
rows = np.array(rows)
a, b = np.polyfit(rows[:, 0], rows[:, 1], 1)
print(f"  fit H(gap) = {a:.3f} log2(ln x) + {b:.3f}  (a random-exponential model with mean ln x gives slope 1)")

print("E-F2 Benford: leading digit of the primes below 1e8 against Benford's log10(1 + 1/d) and against uniform-in-first-digit:")
d = np.array([int(str(p)[0]) for p in pr[::50]])
cnt = np.bincount(d, minlength=10)[1:] / len(d)
benford = np.log10(1 + 1 / np.arange(1, 10))
print("  digit:   " + " ".join(f"{k:6d}" for k in range(1, 10)))
print("  primes:  " + " ".join(f"{v:6.4f}" for v in cnt))
print("  Benford: " + " ".join(f"{v:6.4f}" for v in benford))
print(f"  primes are nearly flat (Luque-Lacasa 2009: generalized Benford with exponent alpha ~ 1/(ln N - 1.1) = {1/(log(X)-1.1):.3f}); chi-square to Benford {((cnt-benford)**2/benford).sum()*len(d):.0f} on {len(d)} primes")

print("E-F5 digit-sum parity of primes (Mauduit-Rivat: equidistributed): excess of even digit sums per decade")
s = np.array([sum(map(int, str(p))) for p in pr[::10]])
for lo, hi in ((10**4, 10**5), (10**5, 10**6), (10**6, 10**7), (10**7, 10**8)):
    sub = pr[::10]; m = (sub >= lo) & (sub < hi)
    ev = (s[m] % 2 == 0).mean()
    n = m.sum()
    print(f"  [{lo:.0e},{hi:.0e}): P(even digit sum) = {ev:.4f} +- {sqrt(0.25/n):.4f} (n = {n} sampled)")
print(f"[{time.time()-t0:.0f} s]")
