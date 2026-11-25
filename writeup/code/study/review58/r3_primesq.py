"""Expected contribution of the prime squares (and fourth powers) to G(T) for a symplectic family (mean sum alpha^2 = -1,
mean sum alpha^4 = -1 for USp(4); sum alpha^2 mean -1 for Sym^3 of SU(2) (USp(4) image), +1 for RS (orthogonal)).
K2(T) = (1/2pi) int_0^T e^{2i theta(t)} * m2 * 2 Re P1(1 + 2it) dt,  P1(s) = sum_p ln p p^-s = -zeta'/zeta(s) - sum_{k>=2} sum_p ln p p^-ks.
(sign: (L'/L)(1/2+it) prime-square part = -sum_p Lambda(p^2) p^{-1-2it} = -m2 * P1(1+2it) in the mean; plus the dual term.)
Hence K2 = -(m2/2pi) int e^{2i theta} 2 Re P1(1+2it) dt. With m2 = -1: +(1/2pi) int e^{2i theta} 2 Re P1 dt."""
import json, math, sys
import numpy as np, mpmath
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
mpmath.mp.dps = 20
dt = 0.005; T = np.arange(dt / 2, 150, dt)
dc = 0.02; Tc = np.arange(0, 150 + 2 * dc, dc)
N = 10 ** 4; s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, 101):
    if s[i]: s[i * i::i] = False
pr = np.nonzero(s)[0].astype(float); lp = np.log(pr)
import os
fn = 'r3_P1.npz'
if os.path.exists(fn):
    Z = np.load(fn); Pa, Pb = Z['Pa'], Z['Pb']
else:
    Pa = np.zeros(len(Tc), complex); Pb = np.zeros(len(Tc), complex)
    def mz(z): return complex(-mpmath.zeta(z, derivative=1) / mpmath.zeta(z))
    for i, t in enumerate(Tc):
        if t == 0:
            Pa[i] = np.nan; Pb[i] = np.nan; continue
        # P1(s) = sum_k mu(k) (-zeta'/zeta)(k s)
        ex = np.exp(-(1 + 2j * t) * lp)
        hp = sum(np.sum(lp * ex ** m) for m in range(2, 15))
        Pa[i] = mz(mpmath.mpc(1, 2 * t)) - hp
        Pb[i] = np.sum(lp * ex ** 2)
    np.savez(fn, Pa=Pa, Pb=Pb)
# check P1 against the direct prime sum at a few points (k>=2 part only is checkable directly; check 2+4it)
for t in (1.0, 7.3, 40.0):
    i = int(round(t / dc)); direct = np.sum(lp * np.exp(-(2 + 4j * Tc[i]) * lp))
    print('check P1(2+4it) t=%.2f stored %s direct(p<=1e4) %s' % (Tc[i], np.round(Pb[i], 6), np.round(direct, 6)))
ReA = np.interp(T, Tc[1:], Pa[1:].real); ReB = np.interp(T, Tc[1:], Pb[1:].real)
P1a = ReA; P1b = ReB
rows = [r for r in load()]
Ts = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
out = {}
for cls, m2, m4 in (('g2', -1, -1), ('sym3', -1, None), ('rs', +1, None)):
    S = [r for r in rows if r['cls'] == cls]
    K = []
    for r in S:
        e = np.exp(2j * theta(T, r['Q'], r['mu']))
        integ = -(m2 / (2 * math.pi)) * e * 2 * P1a * dt
        if m4 is not None: integ = integ - (m4 / (2 * math.pi)) * e * 2 * P1b * dt
        cum = np.cumsum(integ)
        K.append([cum[min(len(T) - 1, int(x / dt) - 1)] for x in Ts])
    K = np.array(K)
    out[cls] = K
    print(cls, 'mean prime-square(+p^4) term Re/Im at T =', Ts)
    print('  Re', np.round(K.real.mean(0), 3))
    print('  Im', np.round(K.imag.mean(0), 3))
np.savez('r3_K.npz', **{k: v for k, v in out.items()})
