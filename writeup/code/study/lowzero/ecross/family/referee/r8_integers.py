"""Referee check 8: Petersson c-decomposition of the INTEGER resonance sum R_1 = sum_n lam_Delta(n) w_n (Hann [20,400]),
whose value -(1/2pi) int H is a theorem (Voronoi). Does c = 1 carry all of it, or only 1/omega as for primes?"""
import math, numpy as np, cypari2
from scipy.special import jv, loggamma, digamma
TP = 2*math.pi
pari = cypari2.Pari(); pari.allocatemem(300*10**6)
X = 5400
tau = np.array([float(int(v)) for v in pari('mfcoefs(mfDelta(),%d)' % X)])
n = np.arange(X+1, dtype=float); lam = np.zeros(X+1); lam[1:] = tau[1:]/n[1:]**5.5
nn = np.arange(2, X+1, dtype=float); t = TP*np.sqrt(nn)
f = lambda t: 2*np.real(digamma(6 + 1j*t)) - 2*math.log(TP) - np.log(nn)
for _ in range(40):
    d = (f(t+1e-5) - f(t-1e-5))/2e-5; t = t - f(t)/d
phi2 = (f(t+1e-5) - f(t-1e-5))/2e-5
phase = 2*(-t*math.log(TP) + np.imag(loggamma(6 + 1j*t))) - t*np.log(nn)
H = np.where((t > 20) & (t < 400), np.sin(math.pi*(t-20)/380)**2, 0.0)
wn = np.zeros(X+1, complex); wn[2:] = -(1/TP)*nn**-0.5*np.sqrt(TP/phi2)*np.exp(1j*(phase + math.pi/4))*H
om = 2.8402873751675
idx = np.arange(2, X+1)
def Srow(c):
    if c == 1: return np.ones(c)
    x = np.arange(1, c, dtype=np.int64); x = x[np.gcd(x, c) == 1]
    inv = np.array([pow(int(a), -1, c) for a in x], dtype=np.int64)
    a = np.arange(c, dtype=np.int64)
    return np.cos(TP*((np.outer(a, x) + inv[None, :]) % c)/c).sum(axis=1)
cum = 0j; R1 = np.sum(lam[idx]*wn[idx]); out = {}
lamrec = np.zeros(len(idx))
for c in range(1, 401):
    term = TP/om*Srow(c)[idx % c]/c*jv(11, 4*math.pi*np.sqrt(idx)/c)
    lamrec += term
    cum += np.sum(term*wn[idx])
    if c in (1, 2, 5, 10, 20, 40, 60, 100, 400): out[c] = cum
print('integers: R_1 = %.3f%+.3fi ; max |lam rebuilt - lam| = %.1e' % (R1.real, R1.imag, np.max(np.abs(lamrec - lam[idx]))))
for c, v in out.items():
    print('  sum_{c<=%3d}: %8.3f%+.3fi  fraction of R_1: %.3f' % (c, v.real, v.imag, v.real/R1.real))
print('  1/omega = %.4f' % (1/om))
