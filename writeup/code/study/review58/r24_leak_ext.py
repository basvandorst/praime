"""Line leakage with lines built from p in {2,3,5,7,11,13} (vectorised); compares Re/Im of the family-mean leakage with the
observed residual and the window over-dispersion after subtraction."""
import json, math, sys
import numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
rows = [r for r in load() if r['cls'] == 'g2']
def cc(x):
    try: return complex(x)
    except Exception: return complex(float(x))
PR = {2: 16, 3: 10, 5: 6, 7: 4, 11: 3, 13: 3}
def Sp(L, p, Vmax):
    P = 1 / pari.lfuneuler(L, p)
    co = [cc(c).real for c in pari.Vecrev(P)] if str(pari.type(P)) == 't_POL' else [cc(P).real]
    mu = [c * p ** (-j / 2) for j, c in enumerate(co)]
    J = Vmax + len(mu) + 2; lam = [1.0] + [0.0] * J
    for j in range(1, J + 1): lam[j] = -sum(mu[i] * lam[j - i] for i in range(1, min(j, len(mu) - 1) + 1))
    F = sum(lam[j] * mu[j] * p ** (-j) for j in range(len(mu)))
    vs = list(range(-(len(mu) - 1), Vmax + 1))
    S = np.array([sum(lam[h + v] * mu[h] * p ** (-(2 * h + v) / 2) for h in range(len(mu)) if h + v >= 0) for v in vs])
    return np.array(vs), S / F
Ts = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150], float)
LE = []; OB = []
for r in rows:
    L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2]))
    lr = np.zeros(1); amp = np.ones(1)
    for p, V in PR.items():
        vs, s = Sp(L, p, V)
        lr = (lr[:, None] + vs[None, :] * math.log(p)).ravel(); amp = (amp[:, None] * s[None, :]).ravel()
    keep = np.abs(lr) > 1e-12; lr = lr[keep]; amp = amp[keep] * (-r['eps'] * r['c'])
    lk = np.array([np.sum(amp / (2 * math.pi) * (np.exp(1j * lr * T) - 1) / (1j * lr)) for T in Ts])
    LE.append(lk)
    g = r['g']; ph = np.exp(2j * theta(g, r['Q'], r['mu']))
    OB.append(np.array([ph[(g > 0) & (g <= T)].sum() + r['eps'] * r['c'] * T / (2 * math.pi) for T in Ts]))
LE = np.array(LE); OB = np.array(OB)
np.savez('r24_leak_ext.npz', LE=LE, OB=OB, Ts=Ts)
for j, T in enumerate(Ts):
    d = OB[:, j] - LE[:, j]
    print('T %3.0f  obs %+5.2f%+5.2fi  leak %+5.2f%+5.2fi  obs-leak Re %+5.2f +- %.2f  Im %+5.2f +- %.2f' % (T, OB[:, j].real.mean(), OB[:, j].imag.mean(),
          LE[:, j].real.mean(), LE[:, j].imag.mean(), d.real.mean(), jk_mean(d.real)[1], d.imag.mean(), jk_mean(d.imag)[1]))
rng = np.random.default_rng(9)
for lab, X in (('observed', OB), ('obs - leak(2..13)', OB - LE)):
    W = np.diff(X[:, 1:], axis=1); n = len(W)
    for part in ('real', 'imag'):
        Y = W.real if part == 'real' else W.imag
        st = lambda Z_: np.sum((Z_.mean(0) / (Z_.std(0, ddof=1) / math.sqrt(n))) ** 2)
        obs = st(Y); null = np.array([st(Y * rng.choice([-1, 1], size=(n, 1))) for _ in range(3000)])
        print('%-18s %s: chi2 %.1f on %d windows, null p = %.4f' % (lab, part, obs, Y.shape[1], np.mean(null >= obs)))
