# A(F) = 4 int_0^inf Fhat cosh(pi xi) dxi + 2 int_0^inf F(t) Omega_inf(t) dt for the reconstructed F (expected 0), and F >= 0 up to t = 150.
import numpy as np, mpmath as mp
d = np.load('/tmp/claude-1000/conv_prover/fhat.npz'); xis = d['xi']; Fh = d['Fh']; h = xis[1]-xis[0]
w = np.full(xis.size, h); w[0] = h/2
pole = 4*np.sum(w*Fh*np.cosh(np.pi*xis))
intF = 2*np.sum(w*Fh)      # = F(0)... no: int F dt = Fhat(0)
dt = 0.01; ts = np.arange(0, 150+dt/2, dt)
Ft = np.array([2*np.sum(w*Fh*np.cos(2*np.pi*xis*t)) for t in ts])
Om = np.array([float(mp.re(mp.digamma(mp.mpf(1)/4 + 1j*mp.mpf(t)/2)) - mp.log(mp.pi))/(2*np.pi) for t in ts])
wt = np.full(ts.size, dt); wt[0] = dt/2
arch = 2*np.sum(wt*Ft*Om)
print("int F dt (trapz) = %.12e   Fhat(0) = %.12e" % (2*np.sum(wt*Ft), Fh[0]))
print("pole part = %.12e  arch part = %.12e  A(F) = %.3e   A(F)/int F = %.3e" % (pole, arch, pole+arch, (pole+arch)/Fh[0]))
print("min F on [0,150] = %.3e ; max|F| on [100,150] = %.3e ; F(0) = %.3e" % (Ft.min(), np.abs(Ft[ts>=100]).max(), Ft[0]))
neg = ts[Ft < -1e-14]
print("points with F < -1e-14:", neg[:10], len(neg))
