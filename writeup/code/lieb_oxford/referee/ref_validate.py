"""Validate the referee's Psi against (i) Monte Carlo of definition (19), (ii) the paper's closed forms
(43) (sphere) and (78) (ball), read from the LaTeX source."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ref_core import W_ss, W_BB, Psi_ball, Phi_sphere

rng = np.random.default_rng(12345)

def unif_ball(n, R):
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    return v * (R * rng.random(n) ** (1 / 3))[:, None]

def unif_sphere(n, R):
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    return v * R

def mc_W(R1, R2, kind, n=4_000_000, reps=5):
    e1 = np.array([1.0, 0, 0])
    vals = []
    for _ in range(reps):
        if kind == 'ball':
            u = unif_ball(n, R1); v = unif_ball(n, R2)
        else:
            u = unif_sphere(n, R1); v = unif_sphere(n, R2)
        d = np.linalg.norm(e1 + u - v, axis=1)
        vals.append(np.mean(1.0 / d))
    vals = np.array(vals)
    return vals.mean(), vals.std(ddof=1) / np.sqrt(reps)

print("--- Monte Carlo of definition (19) vs referee closed/shell forms")
cases = [(0.5, 0.3), (2.0, 0.7), (0.4, 1.5), (3.0, 2.5), (1.2, 0.2), (0.9, 0.9), (0.2, 4.0), (0.6, 0.45)]
for R1, R2 in cases:
    # sphere case: 1/|.| integrable singularity on intersecting spheres -> MC variance large but finite? use ball case mainly
    m, s = mc_W(R1, R2, 'ball', n=2_000_000, reps=5)
    w = float(W_BB(R1, R2))
    print(f"ball  R1={R1:4} R2={R2:4}  MC {m:.6f} +- {s:.1e}   shell {w:.6f}   z = {(w-m)/s:+.2f}")
for R1, R2 in cases:
    m, s = mc_W(R1, R2, 'sphere', n=2_000_000, reps=5)
    w = float(W_ss(R1, R2))
    print(f"sphere R1={R1:4} R2={R2:4}  MC {m:.6f} +- {s:.1e}   closed {w:.6f}   z = {(w-m)/s:+.2f}")

# ---- paper formulas
def psi78(a, b):
    x1 = np.maximum(a + b - a * b, 0.0); x2 = np.maximum(np.abs(a - b) - a * b, 0.0)
    return (x1 ** 4 / 160 * (a*a*b*b - 5*a*a - 5*b*b + 4*a*b*b + 4*b*a*a + 20*a*b)
            - x2 ** 4 / 160 * (a*a*b*b - 5*a*a - 5*b*b + 4*np.abs(a - b)*a*b - 20*a*b))

def psi43(a, b):
    return a*a*b*b/4 * (np.maximum(a + b - a*b, 0)**2 - np.maximum(np.abs(a - b) - a*b, 0)**2)

print("--- referee Psi vs paper (78) and (43) at 200000 random points (relative to a^3 b^3)")
a = np.exp(rng.uniform(np.log(1e-3), np.log(1e3), 200000))
b = np.exp(rng.uniform(np.log(1e-3), np.log(1e1), 200000))
mine = Psi_ball(a, b); paper = psi78(a, b)
err = np.abs(mine - paper) / (a**3 * b**3)
print("ball  : max |mine - (78)| / a^3b^3 =", err.max(), " at", a[err.argmax()], b[err.argmax()])
mine = Phi_sphere(a, b, 1.0, 1.0); paper = psi43(a, b)
err = np.abs(mine - paper) / (a**3 * b**3)
print("sphere: max |mine - (43)| / a^3b^3 =", err.max(), " at", a[err.argmax()], b[err.argmax()])
# special lines
for (aa, bb) in [(3.0, 1.5), (2.0, 2.0), (1.0, 0.5), (0.5, 1.0), (7/8, 7/8), (100.0, 0.618), (0.3, 0.3/1.3)]:
    print(f"  (a,b)=({aa:.4f},{bb:.4f}): mine {float(Psi_ball(aa, bb)):.15e}  (78) {float(psi78(aa, bb)):.15e}")
