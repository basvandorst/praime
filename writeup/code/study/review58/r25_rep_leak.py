"""Replication with the line leakage (p <= 13) subtracted: residual Re/Im(G - pred - leak) for the replication curves (rep/*.npz),
by root number, profile in T, and pooled with rounds 151+152 (r24_leak_ext.npz)."""
import json, math, sys, glob
import numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import theta, jk_mean
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
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
rep = []
for fn in sorted(glob.glob('rep/*.npz')):
    z = np.load(fn); info = json.loads(str(z['info']))
    if info.get('status') != 'COMPLETE': continue
    rep.append((info, z['g']))
LE = []; OB = []; EP = []; CC = []; NN = []
for info, g in rep:
    f, h = info['spec'][1], info['spec'][2]
    L = pari('lfungenus2([%s,%s])' % (f, h)); eps = info['eps'][0]; c = info['c']
    lr = np.zeros(1); amp = np.ones(1)
    for p, V in PR.items():
        vs, s = Sp(L, p, V)
        lr = (lr[:, None] + vs[None, :] * math.log(p)).ravel(); amp = (amp[:, None] * s[None, :]).ravel()
    keep = np.abs(lr) > 1e-12; lr = lr[keep]; amp = amp[keep] * (-eps * c)
    LE.append(np.array([np.sum(amp / (2 * math.pi) * (np.exp(1j * lr * T) - 1) / (1j * lr)) for T in Ts]))
    ph = np.exp(2j * theta(g, info['Q'], info['mu']))
    OB.append(np.array([ph[(g > 0) & (g <= T)].sum() + eps * c * T / (2 * math.pi) for T in Ts]))
    EP.append(eps); CC.append(c); NN.append(info['Q'])
LE = np.array(LE); OB = np.array(OB); EP = np.array(EP); CC = np.array(CC)
np.savez('r25_rep_leak.npz', LE=LE, OB=OB, EP=EP, CC=CC, NN=np.array(NN), Ts=Ts)
Z = np.load('r24_leak_ext.npz'); LE0, OB0 = Z['LE'], Z['OB']
print('replication n = %d (eps +1: %d, mean c %.2f; eps -1: %d, mean c %.2f)' % (len(OB), (EP > 0).sum(), CC[EP > 0].mean(), (EP < 0).sum(), CC[EP < 0].mean()))
for lab, O, Lk in (('replication', OB, LE), ('rounds 151+152', OB0, LE0), ('pooled', np.r_[OB0, OB], np.r_[LE0, LE])):
    D = O - Lk
    print('%-15s ' % lab + ' '.join('T%d %+.2f(%.2f)' % (T, *jk_mean(D[:, j].real)) for j, T in enumerate(Ts) if T in (20, 30, 60, 90, 120, 150)))
    print('%-15s raw Re(T=150) %+.2f +- %.2f ; minus leak %+.2f +- %.2f ; Im minus leak %+.2f +- %.2f ; increment (30,150] minus leak %+.2f +- %.2f' % (
        '', *jk_mean(O[:, -1].real), *jk_mean(D[:, -1].real), *jk_mean(D[:, -1].imag), *jk_mean((D[:, -1] - D[:, 2]).real)))
s = EP > 0
print('replication by eps (minus leak, T=150): +1 %+.2f +- %.2f ; -1 %+.2f +- %.2f' % (*jk_mean((OB - LE)[s, -1].real), *jk_mean((OB - LE)[~s, -1].real)))
two = np.array([int(round(n)) % 2 == 0 for n in NN])
if two.any(): print('replication 2 | N: n=%d mean raw %+.2f ; odd N: %+.2f' % (two.sum(), OB[two, -1].real.mean(), OB[~two, -1].real.mean()))
# registered-type slopes on the replication
P = -EP * CC * 150 / (2 * math.pi); G = OB[:, -1] + P * 0 - (-P) * 0
Gfull = OB[:, -1] + (-EP * CC * 150 / (2 * math.pi))      # residual + pred = G
for nm, sel in (('all', np.ones(len(EP), bool)), ('eps+1', s), ('eps-1', ~s)):
    a = Gfull[sel].real * P[sel]; b = P[sel] ** 2; sx, sy = a.sum(), b.sum(); loo = (sx - a) / (sy - b); n = sel.sum()
    print('replication slope %-6s %.3f +- %.3f' % (nm, sx / sy, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))))
