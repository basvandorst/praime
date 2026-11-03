import sys, time
import numpy as np
from math import log, sqrt, pi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, unfold, spectrum, main_term, C2

t0 = time.time()
X = 10**8
isp = sieve(X + 2)
pr = np.flatnonzero(isp[10**7:X]) + 10**7


def digit_sum_base(n, b):
    s = np.zeros_like(n); m = n.copy()
    while m.any():
        s += m % b; m //= b
    return s


print("E-F6 digit-sum parity in other bases, primes in [1e7, 1e8) against integers coprime to 6 and to 2310 (mechanism: b = -1 mod (b+1), so the alternating digit sum is n mod (b+1)):")
n = np.arange(10**7, 10**8, 3)
for b, label in ((2, "base 2 (b+1 = 3)"), (3, "base 3 (b+1 = 4)"), (4, "base 4 (b+1 = 5)"), (6, "base 6 (b+1 = 7)"), (10, "base 10 (b+1 = 11)")):
    sp = digit_sum_base(pr, b) % 2 == 0
    keep6 = (n % 2 == 1) & (n % 3 != 0)
    keep2310 = keep6 & (n % 5 != 0) & (n % 7 != 0) & (n % 11 != 0)
    s6 = digit_sum_base(n[keep6], b) % 2 == 0
    s2310 = digit_sum_base(n[keep2310], b) % 2 == 0
    print(f"  {label:18s}: primes {sp.mean():.4f} +- {sqrt(0.25/len(pr)):.4f} | coprime to 6 {s6.mean():.4f} | coprime to 2310 {s2310.mean():.4f} +- {sqrt(0.25/keep2310.sum()):.4f}  [{time.time()-t0:.0f} s]", flush=True)

print("E-B9 Lehmer pairs: unfolded spacings below 0.05 per 10000 zeros, against GUE (pi^2/9) s0^3 * 10000 = 1.37 and Poisson 488:")
gam = zeros()
s = unfold(gam)
for k in range(0, 100000, 20000):
    seg = s[k:k + 20000]
    print(f"  zeros {k}..{k+20000}: {(seg < 0.05).sum()} pairs below 0.05 (GUE 2.7 per 20000), {(seg < 0.1).sum()} below 0.1 (GUE 21.9)")

print("E-A10 blind zero finding from the twin primes alone (twins in [1e3, 1e8], weights ln^2 p / 2C2, scan [10, 40] step 0.02):")
L = 10**3
prall = np.flatnonzero(isp[:X + 1]); prall = prall[prall >= L]
tw = prall[isp[prall + 2]].astype(np.float64)
tg = np.arange(10.0, 40.0, 0.02)
u = log(X / L)
S = (spectrum(tw, np.log(tw) ** 2 / (2 * C2), tg) - main_term(X, L, tg)).real / u
found = []
for i in range(1, len(tg) - 1):
    if S[i] < S[i - 1] and S[i] < S[i + 1] and S[i] < -1.0:
        a, b, c = S[i - 1], S[i], S[i + 1]
        found.append(tg[i] + 0.5 * (a - c) / (a - 2 * b + c) * 0.02)
found = np.array(found)
true = gam[(gam > 10) & (gam < 40)]
hit = np.array([np.abs(found - g).min() for g in true]) if len(found) else np.full(len(true), np.inf)
print(f"  {len(found)} dips below -1.0, {len(true)} true zeros; zeros with a dip within 0.1: {(hit < 0.1).sum()}/{len(true)}, spurious dips: {sum(np.abs(true - f).min() > 0.1 for f in found)}; depths at the true zeros: " + ", ".join(f"{S[np.argmin(np.abs(tg - g))]:.2f}" for g in true))
print(f"[{time.time()-t0:.0f} s]")
