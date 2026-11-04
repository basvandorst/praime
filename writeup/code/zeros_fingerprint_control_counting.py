import sys, time
import numpy as np
from math import log, sqrt, pi, e
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, sieve, primes_upto, spectrum, psi_from_zeros, lam_array

t0 = time.time()
gam = zeros()
print("E-B14 control: the fingerprint on a GUE spectrum rescaled to the zeros' density (no arithmetic inside): mean cos(x ln m) for m = 2, 3, 5 against the zeros' -Lambda(m)/(sqrt m ln(T/2pi))")
rng = np.random.default_rng(4)
N = 4000
A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)); H = (A + A.conj().T) / 2
ev = np.sort(np.linalg.eigvalsh(H))[N // 4: 3 * N // 4]
cdf = (np.arange(len(ev)) + 0.5)
fake = np.interp(cdf, (np.arange(1, 100001) - 0.5)[:len(ev)], gam[:len(ev)]) if False else None
g = gam[:len(ev)]
dens = np.log(g / (2 * pi)) / (2 * pi)
u = np.cumsum(np.concatenate([[0.0], np.diff(ev) / np.diff(ev).mean()]))
fake = np.interp(u, np.cumsum(np.concatenate([[0.0], np.diff(g) * dens[:-1]])), g)
for m in (2, 3, 5, 7):
    print(f"  m={m}: GUE-shaped sequence {np.cos(fake * log(m)).mean():+.4f}   zeros (same count) {np.cos(g * log(m)).mean():+.4f}   law {-log(m)/sqrt(m)/log(g.mean()/(2*pi)):+.4f}")

print("E-B15 the fingerprint with only the first 1000 zeros (T = 1419): mean cos(gamma ln m) against -Lambda(m)/(sqrt m ln(T/2pi)):")
g1 = gam[:1000]
print("  " + "  ".join(f"m={m}: {np.cos(g1*log(m)).mean():+.3f} ({-(log(2) if m in (2,4,8) else log(3) if m in (3,9) else (0 if m in (6,10) else log(m)))/sqrt(m)/log(g1.mean()/(2*pi)):+.3f})" for m in (2, 3, 4, 5, 6, 7, 10, 11)))

print("E-A15 blind zeros from the Liouville function alone (lambda(n), n <= 1e6): peaks of |sum lambda(n) n^(-1/2-it)| on [10, 60]")
lam = np.ones(10**6 + 1, dtype=np.int8); lam[0] = 0
for q in primes_upto(10**6):
    qq = q
    while qq <= 10**6:
        lam[qq::qq] *= -1; qq *= q
nn = np.arange(1, 10**6 + 1, dtype=np.float64)
tg = np.arange(10.0, 60.0, 0.02)
S = np.abs(spectrum(nn, lam[1:].astype(np.float64), tg))
found = []
for i in range(1, len(tg) - 1):
    if S[i] > S[i-1] and S[i] > S[i+1] and S[i] > 1.5 * np.median(S):
        a, b, c = S[i-1], S[i], S[i+1]
        found.append(tg[i] + 0.5 * (a - c) / (a - 2 * b + c) * 0.02)
found = np.array(found)
true = gam[(gam > 10) & (gam < 60)]
d = np.array([np.abs(found - z).min() for z in true])
print(f"  {len(found)} peaks, zeros within 0.1: {(d < 0.1).sum()}/{len(true)}, rms {sqrt((d[d<0.5]**2).mean()):.4f}, spurious {sum(np.abs(true - f).min() > 0.1 for f in found)} (zeta(2s)/zeta(s) has poles at the zeros)  [{time.time()-t0:.0f} s]")

print("E-A16 counting primes exactly from the zeros: pi(x) from the 1e5-zero psi (psi -> theta -> pi by rounding the staircase) at x = 1e4, 1e5, 1e6, 1e7:")
isp = sieve(10**7 + 1)
pic = np.cumsum(isp)
for x in (10**4, 10**5, 10**6, 10**7):
    xs = np.array([x + 0.5])
    ps = psi_from_zeros(xs, gam)[0]
    th = ps - sum(psi_from_zeros(np.array([x ** (1 / k) + 0.5]), gam)[0] for k in range(2, 24) if x ** (1 / k) > 2)
    print(f"  x={x:9d}: psi from zeros {ps:14.2f}, theta from zeros {th:14.2f}, true theta {np.log(np.flatnonzero(isp[:x+1]).astype(np.float64)).sum():14.2f}, pi(x) = {pic[x]}  (truncation error ~ x ln x / T = {x*log(x)/gam[-1]:.1f})")
print(f"[{time.time()-t0:.0f} s]")
