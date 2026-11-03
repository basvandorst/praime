import sys, time
import numpy as np
from math import log, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, spectrum, main_term, hl_constant, C2

t0 = time.time()
X = 10**8; L = 10**4
isp = sieve(X + 10)
pr = np.flatnonzero(isp[:X + 1]); pr = pr[pr >= L]
pf = pr.astype(np.float64)
u = log(X / L)
gam = zeros(1001)
z = gam[:1000]; mids = 0.5 * (gam[:300] + gam[1:301])
ts = np.concatenate([z, mids]); main = main_term(X, L, ts)
print("E-D14 triplets and quadruplets at 1000 zeros (p >= 1e4): coefficients k = 3 and 4?")
sub = pf[::7]
rp = (spectrum(sub, 7 * np.log(sub), ts) - main).real / u
dp = rp[:1000].mean() - rp[1000:].mean()
for name, offs, k in (("triplets (0,2,6)", (2, 6), 3), ("triplets (0,4,6)", (4, 6), 3), ("quadruplets (0,2,6,8)", (2, 6, 8), 4)):
    m = np.ones(len(pr), dtype=bool)
    for o in offs: m &= isp[pr + o]
    v = pf[m]; c = hl_constant((0,) + offs)
    r = (spectrum(v, np.log(v) ** k / c, ts) - main).real / u
    d = r[:1000].mean() - r[1000:].mean()
    e = sqrt((r[:1000].std() / sqrt(1000)) ** 2 + (r[1000:].std() / sqrt(300)) ** 2)
    print(f"  {name:22s}: n={len(v):6d}, ratio {d/dp:5.2f} +- {e/abs(dp):.2f} (predicted {k})  [{time.time()-t0:.0f} s]", flush=True)

print("E-Z8 chi3-twisted pairs at the zeros of L(s, chi3): sexy pairs (p, p+6) have chi3(p+6) = chi3(p), so 2; cousins (p, p+4) have p = 1 mod 3 always (chi3 constant), so no L-zero structure (0); (p, p+12): 2")
LZ = np.load("/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins/lzeros.npy", allow_pickle=True)
g3 = LZ[0]; g3 = g3[g3 > 20]
m3 = 0.5 * (g3[:-1] + g3[1:]); ts3 = np.concatenate([g3, m3]); k3 = len(g3)
chi3 = np.where(pr % 3 == 1, 1.0, -1.0)
SP = (spectrum(pf, chi3 * np.log(pf), ts3)).real / u
dP = SP[:k3].mean() - SP[k3:].mean()
print(f"  twisted primes at chi3-zeros {SP[:k3].mean():+.3f}, midpoints {SP[k3:].mean():+.3f}")
for name, h, pred in (("sexy (p, p+6)", 6, 2), ("cousins (p, p+4)", 4, 0), ("(p, p+12)", 12, 2), ("(p, p+10)", 10, "0 (chi3(p+10) = -chi3(p)... = chi3(p+1))")):
    m = isp[pr + h]; v = pf[m]; c = hl_constant([0, h])
    S = (spectrum(v, chi3[m] * np.log(v) * np.log(v + h) / c, ts3)).real / u
    d = S[:k3].mean() - S[k3:].mean()
    e = sqrt((S[:k3].std() / sqrt(k3)) ** 2 + (S[k3:].std() / sqrt(k3 - 1)) ** 2)
    print(f"  twisted {name:16s}: n={len(v):7d}, at L-zeros {S[:k3].mean():+.3f}, midpoints {S[k3:].mean():+.3f}, ratio {d/dP:5.2f} +- {e/abs(dP):.2f} (predicted {pred})  [{time.time()-t0:.0f} s]", flush=True)
