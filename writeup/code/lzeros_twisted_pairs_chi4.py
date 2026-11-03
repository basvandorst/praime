import sys, time
import numpy as np
from math import log, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, spectrum, main_term, hl_constant, C2

t0 = time.time()
LZ = np.load("/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins/lzeros.npy", allow_pickle=True)
g4 = LZ[1]; g4 = g4[g4 > 20]
X = 10**8; L = 10**3
isp = sieve(X + 10)
pr = np.flatnonzero(isp[:X + 1]); pr = pr[pr >= L]
pf = pr.astype(np.float64)
chi4 = np.where(pr % 4 == 1, 1.0, -1.0)
mids = 0.5 * (g4[:-1] + g4[1:])
ts = np.concatenate([g4, mids])
k = len(g4)
u = log(X / L)
print(f"E-Z7 chi4-twisted pair spectra at the {k} zeros of L(s, chi4) above 20 (X = 1e8, p >= 1e3). chi4 is periodic mod 4, so chi4(p+2) = -chi4(p), chi4(p+4) = chi4(p), chi4(p+6) = -chi4(p), chi4(p+8) = chi4(p):")
print("  circle method: the L-zero fluctuation enters once per member with the member's sign, so twisted twins and sexy pairs get 1 - 1 = 0, twisted cousins and (p, p+8) get 1 + 1 = 2, relative to the twisted primes' dip")
SP = (spectrum(pf, chi4 * np.log(pf), ts)).real / u
dP = SP[:k].mean() - SP[k:].mean()
print(f"  twisted primes: at L-zeros {SP[:k].mean():+.3f} +- {SP[:k].std()/sqrt(k):.3f}, midpoints {SP[k:].mean():+.3f}  (predicted -1, 0)  [{time.time()-t0:.0f} s]", flush=True)
for name, h, pred in (("twins", 2, 0), ("cousins", 4, 2), ("sexy", 6, 0), ("(p, p+8)", 8, 2)):
    m = isp[pr + h]
    v = pf[m]
    c = hl_constant([0, h])
    S = (spectrum(v, chi4[m] * np.log(v) * np.log(v + h) / c, ts)).real / u
    d = S[:k].mean() - S[k:].mean()
    e = sqrt((S[:k].std() / sqrt(k)) ** 2 + (S[k:].std() / sqrt(k - 1)) ** 2)
    print(f"  twisted {name:9s}: n={len(v):7d}, at L-zeros {S[:k].mean():+.3f}, midpoints {S[k:].mean():+.3f}, ratio to twisted primes {d/dP:5.2f} +- {e/abs(dP):.2f}  (predicted {pred})  [{time.time()-t0:.0f} s]", flush=True)

print("E-A13 blind recovery of the first zeros of L(s, chi4) and L(s, chi3) from twisted primes alone (primes in [1e3, 1e7], scan [5, 60] step 0.02):")
Xb = 10**7
pb = pr[pr <= Xb].astype(np.float64)
ub = log(Xb / L)
tg = np.arange(5.0, 60.0, 0.02)
for name, chi, true in (("chi4", np.where(pr[pr <= Xb] % 4 == 1, 1.0, -1.0), LZ[1]), ("chi3", np.where(pr[pr <= Xb] % 3 == 1, 1.0, -1.0), LZ[0])):
    S = (spectrum(pb, chi * np.log(pb), tg)).real / ub
    found = []
    for i in range(1, len(tg) - 1):
        if S[i] < S[i - 1] and S[i] < S[i + 1] and S[i] < -0.5:
            a, b, c = S[i - 1], S[i], S[i + 1]
            found.append(tg[i] + 0.5 * (a - c) / (a - 2 * b + c) * 0.02)
    found = np.array(found)
    tr = true[(true > 5) & (true < 60)]
    d = np.array([np.abs(found - z).min() for z in tr])
    print(f"  {name}: {len(found)} dips, {len(tr)} true L-zeros in (5, 60); within 0.1: {(d < 0.1).sum()}, rms {sqrt((d[d < 0.5]**2).mean()):.4f}, spurious {sum(np.abs(tr - f).min() > 0.1 for f in found)}; first four recovered {found[:4].round(4).tolist()} vs {tr[:4].round(4).tolist()}  [{time.time()-t0:.0f} s]")
