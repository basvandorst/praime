"""Second moments M2 = sum^h Z_f^2, M11 = sum^h |Z_f|^2 by Petersson, per modulus c = N m (blocked low-rank kernel).
usage: moment2.py N k T m0:m1 model [lo_fac hi_fac]"""
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, time, numpy as np
from smlib import *
from scipy.special import jv
N, k, T, model = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[5]
M0, M = (int(a) for a in sys.argv[4].split(':'))
lo_fac = float(sys.argv[6]) if len(sys.argv) > 6 else 0.3
hi_fac = float(sys.argv[7]) if len(sys.argv) > 7 else 3.5
out = 'mom_%d_%d_%g_%s_%d_%d_%g.npz' % (N, k, T, model, M0, M, hi_fac)
t0 = time.time()
win = Window(N, k, T, h=0.05)
nlo = int(N*(lo_fac*T/TWO_PI)**2); nhi = int(N*(hi_fac*T/TWO_PI)**2)
n, b = prime_power_coeffs(win, nlo, nhi, model)
ik = (-1)**(k//2)
one = (n == 1)
b1 = complex(b[one][0]) if one.any() else 0j
n2, bb = n[~one], b[~one]
s = TWO_PI*np.sqrt(n2/N)
print('N=%d k=%d T=%g model=%s: %d indices in [%d,%d], beta_1=%.3e%+.3ei; intH=%.5f L=%.5f (%.0fs)' % (N, k, T, model, len(n2), nlo, nhi, b1.real, b1.imag, win.intH, win.L, time.time()-t0), flush=True)
D2 = np.sum(bb*bb); D11 = np.sum(np.abs(bb)**2)
mod = Moduli(N)
G1 = 0j
if M0 == 1 and abs(b1) > 0:
    # G1 = sum_n beta_n Delta_off(1, n) = 2 pi i^{-k} sum_m (1/c) sum_x e(x/c) Ahat(xbar), A(b) = sum_{n = b (c)} beta_n J(4 pi sqrt(n)/c)
    mG = min(int(4*math.pi*math.sqrt(n2.max())/(N*1e-3)) + 1, max(M, 200))
    for m in range(1, mG + 1):
        if m % N == 0 and N >= 1000: continue
        c = N*m
        x, xbar = mod.units_and_bar(m)
        wv = bb*jv(k-1, 4*math.pi*np.sqrt(n2)/c)
        A = np.bincount(n2 % c, weights=wv.real, minlength=c) + 1j*np.bincount(n2 % c, weights=wv.imag, minlength=c)
        Ah = np.fft.ifft(A)*c
        G1 += np.sum(np.exp(2j*math.pi*x/c)*Ah[xbar])/c
    G1 *= TWO_PI*ik
smin, smax = s.min(), s.max()
rows = []
M2c = 0j; M11c = 0j
CH = 5000
ALLM = (N < 1000)          # small levels: include the moduli c = N^2 m' as well
for m in range(M0, M+1):
    if m % N == 0 and not ALLM: continue
    t1 = time.time()
    nc = int(0.6*(smax-smin)*smax/(math.pi*m)) + 40
    nodes = cheb_nodes(smin, smax, nc)
    lam, E = kernel_factor(k, m, nodes, tol=1e-13)
    R = len(lam)
    c = N*m
    x, xbar = mod.units_and_bar(m)
    xm = (-xbar) % c
    res = n2 % c
    arc = np.zeros(len(x), bool); pm1 = (x == 1) | (x == c - 1)
    for q in range(1, m+1):
        if m % q: continue
        d = c//q
        for j in (1, -1):
            arc |= (((x - j) % d) == 0) & (((xbar - j) % d) == 0)
    a2 = 0j; a11 = 0j; p2 = 0j; p11 = 0j
    B = max(1, min(R, int(1.5e8/(16*c))))
    t2 = 0j; t11 = 0j
    for r0 in range(0, R, B):
        r1 = min(R, r0+B)
        A = np.zeros((r1-r0, c), complex)
        for i in range(0, len(s), CH):
            Fc = bary_matrix(nodes, s[i:i+CH]) @ E[:, r0:r1]
            wc = bb[i:i+CH, None]*Fc
            rc = res[i:i+CH]
            for j in range(r1-r0):
                A[j] += np.bincount(rc, weights=wc[:, j].real, minlength=c) + 1j*np.bincount(rc, weights=wc[:, j].imag, minlength=c)
        Ah = np.fft.ifft(A, axis=1)*c
        del A
        for j in range(r1-r0):
            v2 = Ah[j, x]*Ah[j, xbar]; v11 = Ah[j, x]*np.conj(Ah[j, xm])
            t2 += lam[r0+j]*np.sum(v2); t11 += lam[r0+j]*np.sum(v11)
            a2 += lam[r0+j]*np.sum(v2[arc]); a11 += lam[r0+j]*np.sum(v11[arc])
            p2 += lam[r0+j]*np.sum(v2[pm1]); p11 += lam[r0+j]*np.sum(v11[pm1])
        del Ah
    t2 *= TWO_PI*ik/c; t11 *= TWO_PI*ik/c
    a2 *= TWO_PI*ik/c; a11 *= TWO_PI*ik/c; p2 *= TWO_PI*ik/c; p11 *= TWO_PI*ik/c
    M2c += t2; M11c += t11
    rows.append((m, t2.real, t2.imag, t11.real, t11.imag, R, a2.real, a2.imag, a11.real, p2.real, p2.imag, p11.real))
    if m % 10 == 0:
        np.savez(out + '.partial.npz', rows=np.array(rows), D2=D2, D11=D11, b1=b1, G1=G1, L=win.L, intH=win.intH, intH2=win.intH2, N=N, k=k, T=T, M=m)
    if m <= 12 or m % 25 == 0 or m == M:
        print('m=%4d rank %4d: M2 term %+.5f%+.5fi  M11 term %+.5f%+.5fi [arcs %+.4f, x=+-1 %+.4f] | cum M2off %+.5f%+.5fi  M11off %+.5f%+.5fi  (%.1fs)' % (m, R, t2.real, t2.imag, t11.real, t11.imag, a11.real, p11.real, M2c.real, M2c.imag, M11c.real, M11c.imag, time.time()-t1), flush=True)
rows = np.array(rows)
print('diag: D2 = %.5f%+.5fi  D11 = %.5f' % (D2.real, D2.imag, D11), flush=True)
if M0 == 1:
    M2tot = D2 + M2c + b1*b1 + 2*b1*G1
    M11tot = D11 + M11c + abs(b1)**2 + 2*(np.conj(b1)*G1).real
    print('diag: D2 = %.5f%+.5fi  D11 = %.5f ; n=1 terms: b1^2 %.2e, 2 b1 G1 %.2e%+.2ei' % (D2.real, D2.imag, D11, abs(b1)**2, (2*b1*G1).real, (2*b1*G1).imag))
    print('TOTAL (m <= %d): M2 = %.5f%+.5fi   M11 = %.5f%+.5fi ;  L^2 = %.5f, 1.5198 L^2 = %.5f, (1/2pi) int H^2 = %.5f' % (M, M2tot.real, M2tot.imag, M11tot.real, M11tot.imag, win.L**2, 1.5198*win.L**2, win.intH2/TWO_PI))
else:
    print('partial m in [%d,%d]: M2off %+.5f%+.5fi  M11off %+.5f%+.5fi' % (M0, M, M2c.real, M2c.imag, M11c.real, M11c.imag))
np.savez(out, rows=rows, D2=D2, D11=D11, b1=b1, G1=G1, L=win.L, intH=win.intH, intH2=win.intH2, N=N, k=k, T=T, M=M)
print('total %.0fs' % (time.time()-t0))
