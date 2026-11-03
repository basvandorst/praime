import sys, time
import numpy as np
import mpmath as mp
from math import log, sqrt, pi, e
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, sieve, spectrum, main_term, C2, lam_array

mp.mp.dps = 15
TMAX = 260.0
chars = {
    "chi3": (3, [0, 1, -1]),
    "chi4": (4, [0, 1, 0, -1]),
    "chi5r": (5, [0, 1, -1, -1, 1]),
    "chi5c": (5, [0, 1, 1j, -1j, -1]),
}
t0 = time.time()


def lzeros(chi, tmax, step=0.02):
    ts = np.arange(0.5, tmax, step)
    mags = np.array([abs(mp.dirichlet(mp.mpc(0.5, t), chi)) for t in ts], dtype=float)
    out = []
    for i in range(1, len(ts) - 1):
        if mags[i] < mags[i - 1] and mags[i] < mags[i + 1] and mags[i] < 0.25:
            try:
                r = mp.findroot(lambda s: mp.dirichlet(s, chi), mp.mpc(0.5, ts[i]))
            except Exception:
                continue
            if abs(r.real - 0.5) < 1e-6 and abs(float(r.imag) - ts[i]) < 0.1:
                out.append(round(float(r.imag), 7))
    return np.array(sorted(set(out)))


Z = {}
for name, (q, chi) in chars.items():
    Z[name] = lzeros(chi, TMAX)
    g = Z[name]
    expected = (g[-1] / (2 * pi)) * log(q * g[-1] / (2 * pi * e)) + 0.25
    print(f"{name} (mod {q}): {len(g)} zeros below {TMAX} (Riemann-von Mangoldt count ~ {expected:.1f}); first: {', '.join(f'{v:.4f}' for v in g[:5])}  [{time.time()-t0:.0f} s]", flush=True)
np.save("/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins/lzeros.npy", np.array([Z[k] for k in chars], dtype=object), allow_pickle=True)

print("E-Z1 GUE for L-zeros: unfolded spacing variance (density ln(q t / 2pi) / 2pi), GUE 0.178, Poisson 1:")
for name, (q, chi) in chars.items():
    g = Z[name]; g = g[g > 20]
    s = np.diff(g) * np.log(q * g[:-1] / (2 * pi)) / (2 * pi)
    print(f"  {name}: {len(s)} spacings, mean {s.mean():.3f}, var {s.var():.3f}, fraction below 0.25: {(s < 0.25).mean():.3f} (GUE 0.016, Poisson 0.22)")

print("E-Z2 arithmetic fingerprint of L-zeros: mean cos(gamma ln m) predicted -Re(chi(m)) Lambda(m) / (sqrt m ln(q T / 2 pi e)); sign follows the character:")
for name, (q, chi) in chars.items():
    g = Z[name]; T = g[-1]
    norm = log(q * T / (2 * pi * e))
    line = f"  {name}:"
    for m in (2, 3, 4, 5, 7, 8, 9, 11, 13):
        cm = chi[m % q]
        lam = log(2) if m in (2, 4, 8) else (log(3) if m in (3, 9) else log(m))
        pred = -np.real(cm) * lam / (sqrt(m) * norm)
        meas = np.cos(g * log(m)).mean()
        line += f" m={m}: {meas:+.3f} ({pred:+.3f})"
    print(line + f"   [noise ~ {0.7/sqrt(len(g)):.3f}]")

print("E-Z3 superposition: zeta zeros + chi3 zeros below 260, spacing variance of the merged unfolded sequence (independent GUE superposition expected between GUE and Poisson):")
zz = zeros(400); zz = zz[zz < TMAX]
def unf(g, q):
    return np.diff(g) * np.log(q * g[:-1] / (2 * pi)) / (2 * pi)
merged = np.sort(np.concatenate([zz, Z["chi3"]]))
dens = (np.log(merged[:-1] / (2 * pi)) + np.log(3 * merged[:-1] / (2 * pi))) / (2 * pi)
s = np.diff(merged) * dens
print(f"  zeta alone var {unf(zz, 1)[5:].var():.3f}; chi3 alone var {unf(Z['chi3'], 3)[5:].var():.3f}; merged var {s[5:].var():.3f}, fraction below 0.25: {(s[5:] < 0.25).mean():.3f}")
rng = np.random.default_rng(1)
vals = []
for _ in range(3):
    N = 1200
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)); H1 = (A + A.conj().T) / 2
    B = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)); H2 = (B + B.conj().T) / 2
    e1 = np.linalg.eigvalsh(H1)[N//4:3*N//4]; e2 = np.linalg.eigvalsh(H2)[N//4:3*N//4]
    m = np.sort(np.concatenate([e1, e2])); sm = np.diff(m); sm /= sm.mean()
    vals.append((sm.var(), (sm < 0.25).mean()))
vals = np.array(vals)
print(f"  two independent GUE spectra superposed (3 samples): var {vals[:,0].mean():.3f}, fraction below 0.25: {vals[:,1].mean():.3f}")

print("E-Z4 cross-correlation zeta vs chi3 zeros: nearest-distance statistics vs independent expectation:")
d = np.array([np.abs(Z["chi3"] - g).min() for g in zz if g > 20])
dens3 = np.log(3 * zz[zz > 20] / (2 * pi)) / (2 * pi)
print(f"  mean nearest chi3-zero distance from a zeta zero, unfolded by the chi3 density: {(d * dens3).mean():.3f} (independent uniform: 0.25); fraction within 0.05: {(d * dens3 < 0.05).mean():.3f} (independent: 0.10)")

print("E-Z6 twisted spectra at the chi3 zeros (X = 1e8, p >= 1000): primes weighted chi3(p) ln p should dip by -ln(X/L); twins weighted chi3(p) ln^2 p / S(2):")
X = 10**8; L = 10**3
isp = sieve(X + 2)
pr = np.flatnonzero(isp[:X + 1]); pr = pr[pr >= L]
chi = np.array([0, 1, -1])[pr % 3].astype(np.float64)
g3 = Z["chi3"][Z["chi3"] > 20]
mids = 0.5 * (g3[:-1] + g3[1:])
ts = np.concatenate([g3, mids])
u = log(X / L)
pf = pr.astype(np.float64)
SP = spectrum(pf, chi * np.log(pf), ts).real / u
tw = isp[pr + 2]
ST = spectrum(pf[tw], chi[tw] * np.log(pf[tw]) ** 2 / (2 * C2), ts).real / u
k = len(g3)
print(f"  chi3-twisted primes: at L-zeros {SP[:k].mean():+.3f} +- {SP[:k].std()/sqrt(k):.3f}, midpoints {SP[k:].mean():+.3f} +- {SP[k:].std()/sqrt(k-1):.3f}  (predicted -1 at zeros, 0 between)")
print(f"  chi3-twisted twins:  at L-zeros {ST[:k].mean():+.3f} +- {ST[:k].std()/sqrt(k):.3f}, midpoints {ST[k:].mean():+.3f} +- {ST[k:].std()/sqrt(k-1):.3f}  (circle method: coefficient 1 + chi3(2)·(?) ... measured ratio {(ST[:k].mean()-ST[k:].mean())/(SP[:k].mean()-SP[k:].mean()):.2f})")
print(f"[{time.time()-t0:.0f} s]")
