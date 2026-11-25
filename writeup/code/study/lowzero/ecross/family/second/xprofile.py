"""Where in the Kloosterman sum S(n,n';c) = sum_x e((nx + n'xbar)/c) does the m-th off-diagonal term of M11 / M2 live?
Contribution of each unit x, grouped by min(|x|_c, |xbar|_c) (|.|_c = distance to 0 mod c) and by the pairs x = a d +- 1."""
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, time, numpy as np
from smlib import *
N, k, T, model = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
ms = [int(a) for a in sys.argv[5].split(',')]
win = Window(N, k, T, h=0.05)
nlo = int(N*(0.3*T/TWO_PI)**2); nhi = int(N*(3.5*T/TWO_PI)**2)
n, b = prime_power_coeffs(win, nlo, nhi, model)
one = (n == 1); n2, bb = n[~one], b[~one]
s = TWO_PI*np.sqrt(n2/N); smin, smax = s.min(), s.max()
mod = Moduli(N); ik = (-1)**(k//2)
for m in ms:
    nc = int(0.6*(smax-smin)*smax/(math.pi*m)) + 40
    nodes = cheb_nodes(smin, smax, nc); lam, E = kernel_factor(k, m, nodes, tol=1e-13)
    c = N*m; x, xbar = mod.units_and_bar(m); xm = (-xbar) % c; res = n2 % c
    c11 = np.zeros(len(x), complex); c2 = np.zeros(len(x), complex)
    for r in range(len(lam)):
        F = np.concatenate([bary_matrix(nodes, s[i:i+5000]) @ E[:, r] for i in range(0, len(s), 5000)])
        wv = bb*F
        A = np.bincount(res, weights=wv.real, minlength=c) + 1j*np.bincount(res, weights=wv.imag, minlength=c)
        Ah = np.fft.ifft(A)*c
        c11 += lam[r]*Ah[x]*np.conj(Ah[xm]); c2 += lam[r]*Ah[x]*Ah[xbar]
    c11 *= TWO_PI*ik/c; c2 *= TWO_PI*ik/c
    dx = np.minimum(x, c - x); dxb = np.minimum(xbar, c - xbar); sz = np.minimum(dx, dxb)
    # arcs: x = a*(c/q) + j, xbar = a'*(c/q) + j, j = +-1, q | m
    arc = np.zeros(len(x), bool)
    for q in range(1, m+1):
        if m % q: continue
        d = c//q
        for j in (1, -1):
            arc |= (((x - j) % d) == 0) & (((xbar - j) % d) == 0)
    print('m=%d c=%d: M11 term %.5f%+.5fi, M2 term %.5f%+.5fi; #arc units %d' % (m, c, c11.sum().real, c11.sum().imag, c2.sum().real, c2.sum().imag, arc.sum()))
    print('   arcs x = a c/q +-1 (q|m): M11 %.5f%+.5fi  M2 %.5f%+.5fi ; rest: M11 %.5f%+.5fi  M2 %.5f%+.5fi' % (c11[arc].sum().real, c11[arc].sum().imag, c2[arc].sum().real, c2[arc].sum().imag, c11[~arc].sum().real, c11[~arc].sum().imag, c2[~arc].sum().real, c2[~arc].sum().imag))
    for K in [1, 2, 5, 20, 100, 1000]:
        sel = sz <= K
        print('   min(|x|,|xbar|) <= %5d: %6d units, M11 %.5f%+.5fi  M2 %.5f%+.5fi' % (K, sel.sum(), c11[sel].sum().real, c11[sel].sum().imag, c2[sel].sum().real, c2[sel].sum().imag))
    # largest single contributions
    o = np.argsort(-np.abs(c11))[:6]
    print('   top |M11| units: ' + ', '.join('x=%d xbar=%d %.4f%+.4fi' % (x[i], xbar[i], c11[i].real, c11[i].imag) for i in o))
    sabs = np.sort(np.abs(c11))[::-1]
    print('   sum |M11 contributions| = %.4f ; rms per unit %.2e ; sqrt(phi(c))*rms = %.4f' % (np.abs(c11).sum(), np.sqrt(np.mean(np.abs(c11)**2)), math.sqrt(len(x))*np.sqrt(np.mean(np.abs(c11)**2))))
