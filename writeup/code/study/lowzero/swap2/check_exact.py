"""J12(sigma) = int H12(u) K(u) du: Euler product cut at P against the exact value
(PARI lfun for L(1 - iu, psi*) of the primitive character inducing psi, the removed Euler factors at p | q1 q2, p not | f,
and the absolutely convergent A12 to 2e6). No zeros are used. usage: check_exact.py [nreal] [seed]"""
import sys, math, time, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hlib import *

nreal = int(sys.argv[1]) if len(sys.argv) > 1 else 30
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
SIG = (0.5, 1.0)
pari = pari_instance()
t0 = time.time()
PA = 2 * 10 ** 6
prA = primes_upto(PA)
ugs = np.concatenate([-UG[::-1][:-1], UG])            # symmetric grid
wsym = {s: np.exp(-ugs ** 2 / (2 * s * s)) * DU for s in SIG}


def char_spec(c):
    """PARI expression for (G, chi) of a character spec: ('R', D) real Kronecker, ('C', N, m, conj)"""
    if c[0] == 'R':
        return 'znchar(%d)' % c[1]
    N, m, cj = c[1], c[2], c[3]
    x = 'znconreylog(znstar(%d,1),%d)' % (N, m)
    return '[znstar(%d,1), %s]' % (N, ('charconj(znstar(%d,1),%s)' % (N, x)) if cj else x)


def modulus(c):
    return abs(c[1]) if c[0] == 'R' else c[1]


def values(c, pr):
    if c[0] == 'R':
        _, X = real_table(pari, [c[1]], int(pr[-1])); return X[0].astype(complex)
    _, v = conrey_table(pari, c[1], c[2], int(pr[-1]))
    return np.conj(v) if c[3] else v


def exact_J(c1, c2):
    q1, q2 = modulus(c1), modulus(c2)
    pari('A1=%s; A2=%s; NN=lcm(%d,%d); GG=znstar(NN,1)' % (char_spec(c1), char_spec(c2), q1, q2))
    pari('xa=zncharinduce(A1[1],A1[2],NN); xb=zncharinduce(A2[1],A2[2],NN); xp=charmul(GG,xa,charconj(GG,xb)); pp=znchartoprimitive(GG,xp)')
    f = int(pari('pp[1].mod'))
    pari('LL=lfuninit(lfuncreate(pp),[1,0,6.5])')
    Lv = np.array([complex(pari('lfun(LL,1-I*(%.17g))' % u)) for u in ugs])
    H = np.abs(Lv) ** 2
    # removed Euler factors: p | q1 q2, p not | f
    for p in sorted(set(int(x) for x in pari('factor(%d)[,1]' % (q1 * q2)))):
        if f % p == 0: continue
        v = float(pari('chareval(pp[1],pp[2],%d)' % p)); ps = np.exp(2j * math.pi * v)
        H = H * np.abs(1 - ps * np.exp(1j * ugs * math.log(p)) / p) ** 2
    psi = values(c1, prA) * np.conj(values(c2, prA))
    keep = psi != 0; P = prA[keep].astype(float); s = psi[keep]; lp = np.log(P)
    A = np.empty(len(ugs), complex)
    for k, u in enumerate(ugs):
        y = s * np.exp(1j * u * lp) / P
        A[k] = np.exp(np.sum(np.log((1 - 2 / P + 2 * y / P - y * y) / (1 - 1 / P) ** 2)))
    H = H * A
    return {s_: complex(np.sum(wsym[s_] * H)) for s_ in SIG}, f, H


def euler_J(c1, c2, P):
    pr = primes_upto(P)
    psi = values(c1, pr) * np.conj(values(c2, pr))
    return J_pair_direct(psi, pr, SIG)


rng = np.random.default_rng(seed)
# real discriminants of the reserved real designs (keys only) and the r121 complex labels (keys only)
from itertools import combinations
Ds = []
for D in range(-1000, 1001):
    if 100 <= abs(D) <= 1000 and int(pari('isfundamental(%d)' % D)):
        Ds.append(D)
cplx = [(101, 2), (124, 3), (115, 2), (388, 3), (139, 2), (148, 3), (397, 2), (116, 3)]
pairs = []
for _ in range(nreal):
    a, b = rng.choice(len(Ds), 2, replace=False); pairs.append((('R', Ds[a]), ('R', Ds[b])))
pairs += [(('R', -104), ('R', 104)), (('R', 401), ('R', -403)), (('R', -995), ('R', 997))]
for (N1, m1), (N2, m2) in [(cplx[0], cplx[1]), (cplx[2], cplx[3]), (cplx[4], cplx[5]), (cplx[6], cplx[7])]:
    pairs.append((('C', N1, m1, 0), ('C', N2, m2, 0)))
    pairs.append((('C', N1, m1, 0), ('C', N2, m2, 1)))
for N, m in cplx[:4]:
    pairs.append((('C', N, m, 0), ('C', N, m, 1)))           # (chi, conj chi)
for k in range(4):
    pairs.append((('C',) + cplx[k] + (0,), ('R', Ds[rng.integers(len(Ds))])))
rows = []
for c1, c2 in pairs:
    ex, f, H = exact_J(c1, c2)
    row = {'c1': list(c1), 'c2': list(c2), 'f': f}
    for s in SIG: row['exact_%.1f' % s] = [ex[s].real, ex[s].imag]
    for P in (3000, 30000, 300000):
        e = euler_J(c1, c2, P)
        for s in SIG: row['P%d_%.1f' % (P, s)] = [e[s].real, e[s].imag]
    rows.append(row)
    print(c1, c2, 'f', f, ' '.join('s=%.1f exact %.6f%+.6fi P3e3 %.2e P3e4 %.2e P3e5 %.2e' % (
        s, ex[s].real, ex[s].imag, abs(complex(*row['P3000_%.1f' % s]) - ex[s]), abs(complex(*row['P30000_%.1f' % s]) - ex[s]),
        abs(complex(*row['P300000_%.1f' % s]) - ex[s])) for s in SIG), '%.0fs' % (time.time() - t0), flush=True)
for s in SIG:
    ex = np.array([complex(*r['exact_%.1f' % s]) for r in rows])
    for P in (3000, 30000, 300000):
        e = np.array([complex(*r['P%d_%.1f' % (P, s)]) for r in rows])
        d = np.abs(e - ex)
        print('sigma %.1f P %6d: max |dJ| %.2e, rms %.2e, max rel %.2e (|J| range %.3f-%.3f)' % (s, P, d.max(), np.sqrt((d ** 2).mean()), (d / np.abs(ex)).max(), np.abs(ex).min(), np.abs(ex).max()))
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'check_exact_out.json'), 'w'))
print('done %.0fs' % (time.time() - t0))
