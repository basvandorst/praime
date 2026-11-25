"""Second moments M2 = sum^h Z_f^2, M11 = sum^h |Z_f|^2 by Petersson, per modulus c = N m.
usage: moment.py N k T M model [lo_fac hi_fac] [out.npz]"""
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, time, numpy as np
from smlib import *
from scipy.special import jv
N, k, T, model = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[5]
M0, M = (int(a) for a in sys.argv[4].split(':')) if ':' in sys.argv[4] else (1, int(sys.argv[4]))
lo_fac = float(sys.argv[6]) if len(sys.argv) > 6 else 0.3
hi_fac = float(sys.argv[7]) if len(sys.argv) > 7 else 3.5
out = sys.argv[8] if len(sys.argv) > 8 else 'mom_%d_%d_%g_%s_%d_%d_%g.npz' % (N, k, T, model, M0, M, hi_fac)
t0 = time.time()
win = Window(N, k, T, h=0.05)
nlo = int(N*(lo_fac*T/TWO_PI)**2); nhi = int(N*(hi_fac*T/TWO_PI)**2)
n, b = prime_power_coeffs(win, nlo, nhi, model)
ik = (-1)**(k//2)                      # i^{-k}
one = (n == 1)
b1 = complex(b[one][0]) if one.any() else 0j
n2, bb = n[~one], b[~one]
s = TWO_PI*np.sqrt(n2/N)
print('N=%d k=%d T=%g model=%s: %d indices in [%d,%d], beta_1=%.3e%+.3ei; intH=%.5f L=%.5f (%.0fs)' % (N, k, T, model, len(n2), nlo, nhi, b1.real, b1.imag, win.intH, win.L, time.time()-t0), flush=True)
D2 = np.sum(bb*bb); D11 = np.sum(np.abs(bb)**2)
mod = Moduli(N)
# linear forms for the n = 1 cross terms: G1 = sum_n beta_n Delta(1,n) (n>=2), off-diagonal only
KN = np.zeros(N)
xN = np.arange(1, N, dtype=np.int64); iN = mod.invN[xN]
for a0 in range(0, N, 500):
    a = np.arange(a0, min(N, a0+500), dtype=np.int64)
    KN[a] = np.cos(TWO_PI*((xN[None, :] + a[:, None]*iN[None, :]) % N)/N).sum(axis=1)
def Kmod(m):
    if m == 1: return np.ones(1)
    x = np.arange(1, m, dtype=np.int64); x = x[np.gcd(x, m) == 1]; ix = mod.invm(m)[x]
    a = np.arange(m, dtype=np.int64)
    return np.cos(TWO_PI*((x[None, :] + a[:, None]*ix[None, :]) % m)/m).sum(axis=1)
G1 = 0j
mG = int(4*math.pi*math.sqrt(n2.max())/(N*1e-3)) + 1
for m in range(1, mG):
    if m % N == 0: continue
    mb = pow(m, -1, N); Nb = pow(N, -1, m) if m > 1 else 0
    S = KN[(n2*mb*mb) % N]*Kmod(m)[(n2*Nb*Nb) % m]
    G1 += np.sum(bb*S*jv(k-1, 4*math.pi*np.sqrt(n2)/(N*m)))/(N*m)
G1 *= TWO_PI*ik
# kernel nodes
smin, smax = s.min(), s.max()
rows = []
M2c = 0j; M11c = 0j
for m in range(M0, M+1):
    if m % N == 0: continue
    t1 = time.time()
    nc = int(0.6*(smax-smin)*smax/(math.pi*m)) + 40
    nodes = cheb_nodes(smin, smax, nc)
    lam, E = kernel_factor(k, m, nodes, tol=1e-13)
    Fr = np.empty((len(s), len(lam)))
    for i in range(0, len(s), 20000):
        Fr[i:i+20000] = bary_matrix(nodes, s[i:i+20000]) @ E
    c = N*m
    x, xbar = mod.units_and_bar(m)
    xneg = (-x) % c
    res = n2 % c
    t2 = 0j; t11 = 0j
    for r in range(len(lam)):
        wv = bb*Fr[:, r]
        A = np.bincount(res, weights=wv.real, minlength=c) + 1j*np.bincount(res, weights=wv.imag, minlength=c)
        Ah = np.fft.ifft(A)*c                     # sum_b A(b) e(bx/c)
        t2 += lam[r]*np.sum(Ah[x]*Ah[xbar])
        # second factor with conj(beta): Chat(y) = conj(Ahat(-y))
        t11 += lam[r]*np.sum(Ah[x]*np.conj(Ah[(-xbar) % c]))
    t2 *= TWO_PI*ik/c; t11 *= TWO_PI*ik/c
    M2c += t2; M11c += t11
    rows.append((m, t2.real, t2.imag, t11.real, t11.imag, len(lam)))
    if m <= 12 or m % 25 == 0 or m == M:
        print('m=%4d rank %3d: M2 term %+.5f%+.5fi  M11 term %+.5f%+.5fi | cum M2off %+.5f%+.5fi  M11off %+.5f%+.5fi  (%.1fs)' % (m, len(lam), t2.real, t2.imag, t11.real, t11.imag, M2c.real, M2c.imag, M11c.real, M11c.imag, time.time()-t1), flush=True)
rows = np.array(rows)
M2 = D2 + M2c + 2*b1*(np.sum(bb*0) + 0) + 0
# n = 1 contributions: M2: b1^2*Delta(1,1) + 2 b1 * sum_n beta_n Delta(1,n); M11: |b1|^2 + 2 Re(conj(b1) G1)
M2tot = D2 + M2c + b1*b1 + 2*b1*G1
M11tot = D11 + M11c + abs(b1)**2 + 2*(np.conj(b1)*G1).real
print('diag: D2 = %.5f%+.5fi  D11 = %.5f ; n=1 terms: b1^2 %.2e, 2 b1 G1 %.2e%+.2ei (G1 = %.3e%+.3ei)' % (D2.real, D2.imag, D11, abs(b1)**2, (2*b1*G1).real, (2*b1*G1).imag, G1.real, G1.imag))
print('TOTAL (m <= %d): M2 = %.5f%+.5fi   M11 = %.5f%+.5fi ;  L^2 = %.5f, 1.5198 L^2 = %.5f, (1/2pi) int H^2 = %.5f' % (M, M2tot.real, M2tot.imag, M11tot.real, M11tot.imag, win.L**2, 1.5198*win.L**2, win.intH2/TWO_PI))
np.savez(out, rows=rows, D2=D2, D11=D11, b1=b1, G1=G1, L=win.L, intH=win.intH, intH2=win.intH2, N=N, k=k, T=T, M=M)
print('total %.0fs' % (time.time()-t0))
