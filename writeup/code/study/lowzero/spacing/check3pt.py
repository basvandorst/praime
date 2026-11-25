"""Task 1(c), three points: does the wave closure reproduce the exact (ratios conjecture, tplib) change of the
three-point function when a prime is added?  R3 = 1 + T2(12) + T2(13) + T2(23) + T3, T2 = -sinc^2 + C (bklib), T3 from
tplib (T3_gue + Q + P).  Closure: R3_with(0,t,s) = E_phi[ N'(0) N'(t) N'(s) R3_base(N(0), N(t), N(s)) ]."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib')
import numpy as np
import loadlibs, closure
bk = loadlibs.load_bklib(); tp = loadlibs.load_tplib(); G = tp.G
PS = tp.PrimeSums(PB=10 ** 4)

def T2(r, L, S):
    r = np.abs(np.asarray(r, float))
    return -np.sinc(r) ** 2 + bk.C_of_r(G, L, tuple(S), r=r)

def R3(u, v, w, L, S, sine=False):
    if sine:
        a, b, c = np.sinc(v - u), np.sinc(w - v), np.sinc(w - u)
        return 1 - a * a - b * b - c * c + 2 * a * b * c
    d = tp.T3_parts(u, v, w, L, tuple(S), PS=PS)
    T3 = tp.T3_gue(u, v, w) + d['Q'] + d['P']
    return 1 + T2(v - u, L, S) + T2(w - u, L, S) + T2(w - v, L, S) + T3

def R2(r, L, S, sine=False):
    return 1 - np.sinc(r) ** 2 if sine else 1 + T2(r, L, S)

confs = [(t, s) for s in (1.0, 1.6, 2.2, 3.0) for t in (0.35, 0.5 * s, s - 0.35) if t < s]
ph = closure.phases(48)
for L in (9.4, 10.7):
    for p, rest in ((2, ()), (3, ()), (5, ()), (2, (3,))):
        Sb = tuple(sorted(rest + (p,))); w = closure.Wave(p, L)
        print('\nL=%.1f  p=%d made good (base S=%s -> %s)' % (L, p, Sb, rest))
        print('   (t, s)      dR3 exact  closure  c-sine | dT3 exact  closure  c-sine | pair-only dR3')
        sq = []
        for (t, s) in confs:
            ex = R3(np.array([0.0]), np.array([t]), np.array([s]), L, rest)[0] - R3(np.array([0.0]), np.array([t]), np.array([s]), L, Sb)[0]
            x = np.array([0.0, t, s]); res = {}
            for tag, sine in (('ex', False), ('sine', True)):
                N0, Nt, Ns = (w.N(xx, ph) for xx in x); J = w.rho(0.0, ph) * w.rho(t, ph) * w.rho(s, ph)
                r3w = np.mean(J * R3(N0, Nt, Ns, L, Sb, sine)); r3b = R3(np.array([0.0]), np.array([t]), np.array([s]), L, Sb, sine)[0]
                # the pair functions of the closure
                pr = []
                for (a, b) in ((0.0, t), (0.0, s), (t, s)):
                    Jab = w.rho(a, ph) * w.rho(b, ph)
                    pr.append(np.mean(Jab * R2(np.abs(w.N(b, ph) - w.N(a, ph)), L, Sb, sine)) - R2(np.array([b - a]), L, Sb, sine)[0])
                res[tag] = (r3w - r3b, (r3w - r3b) - sum(pr), sum(pr))
            pex = sum((R2(np.array([b - a]), L, rest)[0] - R2(np.array([b - a]), L, Sb)[0]) for (a, b) in ((0.0, t), (0.0, s), (t, s)))
            dT3ex = ex - pex
            sq.append((ex, res['ex'][0], res['sine'][0], dT3ex, res['ex'][1], res['sine'][1]))
            print('  (%.2f,%.2f)  %+.5f  %+.5f  %+.5f | %+.5f  %+.5f  %+.5f | %+.5f' % (t, s, ex, res['ex'][0], res['sine'][0], dT3ex, res['ex'][1], res['sine'][1], pex))
        a = np.array(sq)
        for j, nm in ((1, 'closure, exact base'), (2, 'closure, sine base')):
            print('   dR3: rms(%s - exact) = %.5f  (rms exact %.5f);  dT3: rms diff %.5f (rms exact dT3 %.5f)'
                  % (nm, np.sqrt(np.mean((a[:, j] - a[:, 0]) ** 2)), np.sqrt(np.mean(a[:, 0] ** 2)), np.sqrt(np.mean((a[:, j + 3] - a[:, 3]) ** 2)), np.sqrt(np.mean(a[:, 3] ** 2))))
