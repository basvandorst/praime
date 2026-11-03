import sys, time, lzma
import numpy as np
import mpmath as mp
from math import log, sqrt, pi, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, unfold, smooth_N, GAMMA

t0 = time.time()
gam = zeros()
mp.mp.dps = 12
print("E-X5 RH verified by hand to T = 3000: count the sign changes of Hardy's Z(t) on a grid of step 0.1 and compare with N(T) from the zero table (sign changes = zeros on the line; equality with the total count means every zero is on the line)")
T = 3000
ts = np.arange(10.0, T, 0.1)
Z = np.array([float(mp.siegelz(t)) for t in ts])
changes = int((np.sign(Z[1:]) != np.sign(Z[:-1])).sum())
nT = int((gam < T).sum() - (gam < 10).sum())
print(f"  sign changes of Z in (10, {T}): {changes}; zeros in the table in the same range: {nT}; Riemann-von Mangoldt N({T}) = {smooth_N(T):.1f}  [{time.time()-t0:.0f} s]", flush=True)
miss = np.flatnonzero(np.diff(np.searchsorted(gam, ts)) >= 2)
print(f"  grid cells with two table zeros (where a coarse grid can miss a pair of sign changes): {len(miss)}; refined check of those cells with step 0.002 adds {sum(int((np.sign([float(mp.siegelz(u)) for u in np.arange(ts[i], ts[i+1], 0.002)][1:]) != np.sign([float(mp.siegelz(u)) for u in np.arange(ts[i], ts[i+1], 0.002)][:-1])).sum()) - int(np.sign(Z[i]) != np.sign(Z[i+1])) for i in miss)} sign changes")

print("E-X6 Selberg's central limit theorem: ln|zeta(1/2+it)| for t in [1e4, 1e4+500] step 0.25 (2000 samples), mean and variance against (1/2) ln ln T:")
tt = np.arange(10**4, 10**4 + 500, 0.25)
vals = np.array([float(mp.log(abs(mp.zeta(mp.mpc(0.5, t))))) for t in tt])
vals = vals[np.isfinite(vals)]
print(f"  mean {vals.mean():+.3f} (0 expected), var {vals.var():.3f} against (1/2) ln ln T = {0.5*log(log(10**4)):.3f}; skew {((vals-vals.mean())**3).mean()/vals.std()**3:+.3f}, kurt {((vals-vals.mean())**4).mean()/vals.var()**2:.2f}  [{time.time()-t0:.0f} s]", flush=True)

print("E-X7 Keating-Snaith moments of |zeta(1/2+it)| over t in [1e4, 1e4+500]: <|zeta|^2> against ln(T/2pi) + 2gamma - 1... (Ingham: (1/T) int |zeta|^2 = ln(T/2pi) + 2 gamma - 1), <|zeta|^4> against (1/(2 pi^2)) ln^4 T:")
az = np.exp(vals)
print(f"  <|zeta|^2> = {(az**2).mean():.3f}   (ln(T/2pi) + 2gamma - 1 = {log(10**4/(2*pi)) + 2*GAMMA - 1:.3f});   <|zeta|^4> = {(az**4).mean():.1f}   ((1/2pi^2) ln^4 T = {log(10**4)**4/(2*pi**2):.1f})")

print("E-H1 predictability of zero spacings: regress s_{n+1} on s_n, s_{n-1} for the zeros 20000..100000 and for a GUE sample; R^2:")
s = unfold(gam[20000:])
def r2(seq):
    y = seq[2:]; A = np.vstack([seq[1:-1], seq[:-2], np.ones(len(y))]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    return 1 - ((y - A @ coef) ** 2).sum() / ((y - y.mean()) ** 2).sum(), coef
rz, cz = r2(s)
rng = np.random.default_rng(2)
N = 3000
A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)); H = (A + A.conj().T) / 2
e = np.linalg.eigvalsh(H)[N//4:3*N//4]
dens = np.sqrt(4 * N - 2 * e * e) / (2 * pi) * sqrt(2)
sg = np.diff(e) * dens[:-1]; sg /= sg.mean()
rg, cg = r2(sg)
print(f"  zeros: R^2 {rz:.4f}, coefficients {cz[:2].round(3).tolist()};  GUE sample: R^2 {rg:.4f}, coefficients {cg[:2].round(3).tolist()}")

print("E-H3 compressibility of the spacing sequence (quantised to 32 levels) against its entropy:")
q = np.minimum((s / 3.0 * 32).astype(int), 31).astype(np.uint8)
cnt = np.bincount(q, minlength=32) / len(q)
Hq = -(cnt[cnt > 0] * np.log2(cnt[cnt > 0])).sum()
raw = np.packbits(np.unpackbits(q.reshape(-1, 1), axis=1)[:, 3:].flatten()).tobytes()
comp = lzma.compress(q.tobytes(), preset=9 | lzma.PRESET_EXTREME)
print(f"  entropy {Hq:.3f} bits/spacing = {Hq*len(q)/8/1024:.0f} KB; xz on the byte sequence {len(comp)/1024:.0f} KB; conditional entropy H(s_n | s_(n-1)) = {(-(np.bincount(q[:-1].astype(int)*32+q[1:].astype(int), minlength=1024)/len(q[1:]) * np.log2(np.maximum(np.bincount(q[:-1].astype(int)*32+q[1:].astype(int), minlength=1024)/len(q[1:]), 1e-12))).sum() + (cnt[cnt>0]*np.log2(cnt[cnt>0])).sum()):.3f} bits")

print(f"[{time.time()-t0:.0f} s]")
