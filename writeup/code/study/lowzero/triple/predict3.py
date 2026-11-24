"""Frozen per-L-function predictions of kappa3(l) (third cumulant of the window count, l = 0.5 .. 8) for the zeros of
L(s, chi_D), averaged over the height band with the zero-density weight (4 Gauss nodes in t), under:
  full : ratios conjecture / Bogomolny-Keating, all arithmetic terms, primes | D deleted (tplib)
  gue  : kappa3 of the sine kernel
  Q0   : full, but the bad primes are NOT removed from the diagonal triple-orbit term Q (k[3]_Q over all primes)
  twopt: the connected three-point function of k = 1 at the same L (bad primes act through Sigma2 only)
  base : kappa3 of k = 1 at the same L (no bad prime at all)
  diagonly: base + the diagonal (triple-orbit) change only (bad primes act on kappa3 only through Q)
Also Sigma2 (BK) and k[3] (full) per L-function.
Sets: Fn (fresh: all fundamental D, 1000 < |D| <= 1300, heights 100-600), Bp (B's D at 600-1000, fresh),
      R120, B, Dn (100-600) and C (R120's D at 600-1000): zeros on disk, kappa3 never computed.
usage: predict3.py tables | predict3.py SET (one of Fn Bp R120 B Dn C -> pred3_SET.json) | predict3.py merge -> predictions_kappa3.json"""
import sys, os, json, math, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, 'lib'))
import tplib, bklib
LOW = __import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..', '..', '..')) + '/writeup/code/study/lowzero'
P134 = json.load(open(LOW + '/bk/predictions_per_lfunction.json')); P136 = json.load(open(LOW + '/amp/predictions_alt.json'))
LS = [0.5 * i for i in range(1, 17)]
q = tplib.Quad(); IX = [int(round(l / q.WS)) - 1 for l in LS]
PRIMES = tplib.primes_upto(10 ** 6)


def squarefree(n):
    d = 2
    while d * d <= n:
        if n % (d * d) == 0: return False
        d += 1
    return True


def fundamental(D):
    """D = 1 mod 4 squarefree, or D = 4m with m = 2, 3 mod 4 squarefree (Python's % is non-negative for negative D)"""
    if D % 4 == 1: return squarefree(abs(D))
    if D % 4 == 0:
        m = D // 4
        return m % 4 in (2, 3) and squarefree(abs(m))
    return False


def prime_factors(n):
    n = abs(n); ps = []; p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps


FN = sorted([D for D in range(-1300, 1301) if 1000 < abs(D) <= 1300 and fundamental(D)], key=lambda d: (abs(d), d))
SETS = {'Fn': (FN, 100, 600), 'Bp': (P134['B']['D'], 600, 1000), 'R120': (P134['C']['D'], 100, 600),
        'B': (P134['B']['D'], 100, 600), 'Dn': (P136['Dn']['D'], 100, 600), 'C': (P134['C']['D'], 600, 1000)}

# ---------------- L-grid tables: B1 over all primes at the nodes, Q-part closed form over all primes
t0 = time.time()
Lnodes_all = []
for name, (Ds, lo, hi) in SETS.items():
    for D in Ds:
        Ln, wn = bklib.L_nodes(abs(D), lo, hi, n=4); Lnodes_all += list(Ln)
Lmin, Lmax = min(Lnodes_all) - 0.05, max(Lnodes_all) + 0.05
LGR = np.arange(math.floor(Lmin * 10) / 10, Lmax + 0.1, 0.1)
cache = os.path.join(HERE, 'b1grid.npz'); cdir = os.path.join(HERE, 'b1grid'); os.makedirs(cdir, exist_ok=True)
if os.path.exists(cache) and np.allclose(np.load(cache)['LGR'], LGR):
    z = np.load(cache); B1G = z['B1G']; QG = z['QG']
else:
    PS = tplib.PrimeSums(10 ** 4); B1G = []; QG = []
    for L in LGR:
        f = os.path.join(cdir, 'L%.2f.npz' % L)
        if os.path.exists(f):
            z = np.load(f); B1G.append(z['B1']); QG.append(z['Q']); continue
        if time.time() - t0 > 480: print('time budget used; rerun to continue'); sys.exit(0)
        k = 2 * math.pi / L; arrs = []
        u = 0 * q.T; v = q.T; w = q.S; pos = {'u': u * k, 'v': v * k, 'w': w * k}
        for (A, B1n, B2n) in tplib.TERMS:
            al = -pos[A]; e1 = al + pos[B1n]; e2 = al + pos[B2n]
            arrs.append(PS.B1_all(e1, e2)); arrs.append(PS.B1_all(e2, e1))
        Qv = tplib.kf3_Q_closed(q.lgrid, L, PRIMES); np.savez(f, B1=np.array(arrs), Q=Qv)
        B1G.append(np.array(arrs)); QG.append(Qv)
    B1G = np.array(B1G); QG = np.array(QG); np.savez(cache, LGR=LGR, B1G=B1G, QG=QG)
print('L grid %.1f-%.1f (%d values), B1 tables in %.0fs' % (LGR[0], LGR[-1], len(LGR), time.time() - t0))
MODE = sys.argv[1]
if MODE == 'tables': sys.exit(0)
def term_args(L):
    """the six (f1, f2) argument pairs of B1 at the nodes, in the order used by tplib.T3_parts"""
    k = 2 * math.pi / L; pos = {'u': 0 * q.T * k, 'v': q.T * k, 'w': q.S * k}; out = []
    for (A, B1n, B2n) in tplib.TERMS:
        al = -pos[A]; e1 = al + pos[B1n]; e2 = al + pos[B2n]; out += [(e1, e2), (e2, e1)]
    return out


# B1 has a pole where A_3 has its zero (1 - 2/3 + 3^{-1-x} = 0, x = i pi/ln 3): the p = 3 summand is taken out of the
# interpolated table and added back exactly at each L, so that the pole meets the zero of A_3 in P1 exactly.
B1G3 = B1G.copy()
for i, L in enumerate(LGR):
    for j, (f1, f2) in enumerate(term_args(L)): B1G3[i, j] -= tplib.b1_p(3, f1, f2)
B1S = [[CubicSpline(LGR, B1G3[:, j, :].real, axis=0), CubicSpline(LGR, B1G3[:, j, :].imag, axis=0)] for j in range(6)]
QS = CubicSpline(LGR, QG, axis=0)
kg = q.kf3(tplib.T3_gue(0 * q.T, q.T, q.S)); s2g = tplib.sigma2_gue(q.lgrid); k3g = kg + 3 * s2g - 2 * q.lgrid


def kf3_parts(L, S):
    """k[3] - GUE: (Q-part, P-part) on q.lgrid at fixed L for bad set S."""
    B1all = [B1S[j][0](L) + 1j * B1S[j][1](L) + tplib.b1_p(3, f1, f2) for j, (f1, f2) in enumerate(term_args(L))]
    Qp = QS(L) - (tplib.kf3_Q_closed(q.lgrid, L, S) if S else 0)
    # P-part with the precomputed all-prime B1; Q is handled in closed form, so pass Qall = zeros
    d = tplib.T3_parts(0 * q.T, q.T, q.S, L, S, B1all=B1all, Qall=[np.zeros(len(q.T))] * 3)
    return Qp, q.kf3(d['P'])


def predict(D, lo, hi):
    S = tuple(prime_factors(D)); Ln, wn = bklib.L_nodes(abs(D), lo, hi, n=4)
    acc = {k: 0 for k in ('full', 'Q0', 'twopt', 's2', 'kf3', 'kQ', 'base', 'diagonly')}
    for L, w in zip(Ln, wn):
        Qp, Pp = kf3_parts(L, S); s2 = tplib.sigma2_bk(L, S, q.lgrid)
        kf = kg + Qp + Pp; full = kf + 3 * s2 - 2 * q.lgrid
        Qall = QS(L); Q0 = full + (Qall - Qp)
        if S:
            Qp1, Pp1 = kf3_parts(L, ()); s21 = tplib.sigma2_bk(L, (), q.lgrid)
            twopt = (kg + Qp1 + Pp1) + 3 * s2 - 2 * q.lgrid
            base = (kg + Qp1 + Pp1) + 3 * s21 - 2 * q.lgrid          # kappa3 of k = 1 at the same L
            diagonly = base + (Qp - Qp1)                              # bad primes act only through the triple-orbit term
        else:
            twopt = full; base = full; diagonly = full
        for k_, v in (('full', full), ('Q0', Q0), ('twopt', twopt), ('s2', s2), ('kf3', kf), ('kQ', Qp), ('base', base), ('diagonly', diagonly)):
            acc[k_] = acc[k_] + w * v
    return {k_: [float(x) for x in np.asarray(v)[IX]] for k_, v in acc.items()}, float(np.sum(wn * Ln))


out = {'l': LS, 'gue_kappa3': [float(x) for x in k3g[IX]], 'gue_kf3': [float(x) for x in kg[IX]], 'gue_sigma2': [float(x) for x in s2g[IX]],
       'statistic': 'per L-function: u = (g/2pi) ln(|D| g/(2 pi e)) for the zeros g in [lo, hi]; windows [x, x + l), x from u[0] to u[-1] - l step 0.01; '
                    'c = counts; kappa3 = mean((c - mean c)^3); k[3] = kappa3 - 3 var(c) + 2 mean(c)',
       'model': 'tplib.py (ratios conjecture, Conrey-Snaith 2008 triple correlation, p | k deleted; = Bogomolny-Keating 2013 form to 1e-7), '
                'B1 over p <= 1e4, Q in closed form over p <= 1e6, Sigma2 from the repo bklib (frozen), height average with 4 Gauss nodes'}
t0 = time.time()
if MODE in SETS:
    name = MODE; Ds, lo, hi = SETS[name]
    rows = {'height': [lo, hi], 'D': [int(D) for D in Ds], 'meanL': [], 'full': [], 'Q0': [], 'twopt': [], 'sigma2': [], 'kf3': [], 'kQ': [], 'base': [], 'diagonly': []}
    for D in Ds:
        pr, mL = predict(D, lo, hi)
        rows['meanL'].append(mL); rows['full'].append(pr['full']); rows['Q0'].append(pr['Q0']); rows['twopt'].append(pr['twopt'])
        rows['sigma2'].append(pr['s2']); rows['kf3'].append(pr['kf3']); rows['kQ'].append(pr['kQ']); rows['base'].append(pr['base']); rows['diagonly'].append(pr['diagonly'])
    json.dump(rows, open(os.path.join(HERE, 'pred3_%s.json' % name), 'w'))
    F = np.array(rows['full']); print('%-5s %3d L-functions, heights %d-%d, mean L %.3f (%.0fs); mean full kappa3: %s' % (
        name, len(Ds), lo, hi, np.mean(rows['meanL']), time.time() - t0, ' '.join('%+.4f' % x for x in F.mean(0))))
elif MODE == 'merge':
    for name in SETS: out[name] = json.load(open(os.path.join(HERE, 'pred3_%s.json' % name)))
    txt = json.dumps(out, sort_keys=True)
    open(os.path.join(HERE, 'predictions_kappa3.json'), 'w').write(txt)
    print('sha256 predictions_kappa3.json', hashlib.sha256(txt.encode()).hexdigest())
    print('GUE kappa3 ' + ' '.join('%+.5f' % x for x in k3g[IX]))
