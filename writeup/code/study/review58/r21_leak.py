"""Leakage of the law's own lines into the sharp plain sum.  For a sharp window (0, T], the recipe density of e^{2i theta} has,
besides -eps c/2pi, the off-diagonal terms (A_r/2pi) r^{it}, A_r = -eps c prod_p S_p(v_p)/F_p(1/p) (r = prod p^{v_p}),
S_p(v) = sum_{m-h=v} lambda(p^m) mu(p^h) p^{-(m+h)/2} (real L).  Each contributes (A_r/2pi)(r^{iT} - 1)/(i ln r) to G(T):
bounded in T, but coherent over the family when eps c S_p(v)/F_p has a common sign (it does at p = 2, 3: see r14).
Lines from p in {2,3,5}, -4 <= v_p <= V (v_p < -deg is impossible: S_p(v) = 0).  Writes per-curve leak(T) to r21_leak.npz."""
import json, math, sys
import numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
rows = [r for r in load() if r['cls'] == 'g2']
def cc(x):
    try: return complex(x)
    except Exception: return complex(float(x))
PR = (2, 3, 5); V = {2: 24, 3: 14, 5: 9}
def Sp(L, p, Vmax):
    P = 1 / pari.lfuneuler(L, p)
    co = [cc(c).real for c in pari.Vecrev(P)] if str(pari.type(P)) == 't_POL' else [cc(P).real]
    mu = [c * p ** (-j / 2) for j, c in enumerate(co)]                     # analytic normalisation, k = 2
    J = Vmax + len(mu) + 2; lam = [1.0] + [0.0] * J
    for j in range(1, J + 1): lam[j] = -sum(mu[i] * lam[j - i] for i in range(1, min(j, len(mu) - 1) + 1))
    F = sum(lam[j] * mu[j] * p ** (-j) for j in range(len(mu)))
    S = {}
    for v in range(-(len(mu) - 1), Vmax + 1):
        S[v] = sum(lam[h + v] * mu[h] * p ** (-(2 * h + v) / 2) for h in range(len(mu)) if h + v >= 0)
    return F, S
Ts = np.array([5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150], float)
LEAK = []; OBS = []; info = []
for r in rows:
    L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2]))
    loc = {p: Sp(L, p, V[p]) for p in PR}
    # check: S_p(0) = F_p(1/p)  (the r = 1 local term)
    for p in PR: assert abs(loc[p][1][0] - loc[p][0]) < 1e-9, (r['key'], p)
    ec = r['eps'] * r['c']
    lk = np.zeros(len(Ts), complex); top = []
    for v2 in loc[2][1]:
        for v3 in loc[3][1]:
            for v5 in loc[5][1]:
                if v2 == 0 and v3 == 0 and v5 == 0: continue
                lr = v2 * math.log(2) + v3 * math.log(3) + v5 * math.log(5)
                if abs(lr) < 1e-12: continue
                A = -ec * (loc[2][1][v2] / loc[2][0]) * (loc[3][1][v3] / loc[3][0]) * (loc[5][1][v5] / loc[5][0])
                lk += A / (2 * math.pi) * (np.exp(1j * lr * Ts) - 1) / (1j * lr)
    LEAK.append(lk)
    g = r['g']; ph = np.exp(2j * theta(g, r['Q'], r['mu']))
    OBS.append(np.array([ph[(g > 0) & (g <= T)].sum() + ec * T / (2 * math.pi) for T in Ts]))
    info.append(r['key'])
LEAK = np.array(LEAK); OBS = np.array(OBS)
np.savez('r21_leak.npz', LEAK=LEAK, OBS=OBS, Ts=Ts, keys=np.array(info), eps=np.array([r['eps'] for r in rows]), c=np.array([r['c'] for r in rows]))
print('T     obs Re (jk)      leak Re mean (sd)    obs-leak Re (jk)  |  obs Im   leak Im   obs-leak Im (jk)')
for j, T in enumerate(Ts):
    o = OBS[:, j]; l = LEAK[:, j]; d = o - l
    print('%4.0f  %+5.2f (%.2f)   %+5.2f (%.2f)      %+5.2f (%.2f)    |  %+5.2f   %+5.2f    %+5.2f (%.2f)' % (
        T, o.real.mean(), jk_mean(o.real)[1], l.real.mean(), l.real.std(), d.real.mean(), jk_mean(d.real)[1], o.imag.mean(), l.imag.mean(), d.imag.mean(), jk_mean(d.imag)[1]))
j = len(Ts) - 1
o = OBS[:, j]; l = LEAK[:, j]
print('per-curve corr(obs, leak) at T=150: Re %.3f Im %.3f ; var(obs) %.1f var(obs - leak) %.1f (Re), %.1f %.1f (Im)' % (
    np.corrcoef(o.real, l.real)[0, 1], np.corrcoef(o.imag, l.imag)[0, 1], o.real.var(), (o - l).real.var(), o.imag.var(), (o - l).imag.var()))
for e in (1, -1):
    s = np.array([r['eps'] == e for r in rows])
    print('eps %+d: obs Re %+5.2f  leak Re %+5.2f  obs-leak %+5.2f +- %.2f' % (e, o.real[s].mean(), l.real[s].mean(), (o - l).real[s].mean(), jk_mean((o - l).real[s])[1]))
