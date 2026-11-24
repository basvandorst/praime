"""Data-quality check (not a test statistic): largest unfolded gaps and the count drift N(t) - theta/pi for sets whose
kappa3 residual is an outlier."""
import sys, json, math
import numpy as np
from scipy.special import loggamma
LZ = __import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..', '..', '..')) + '/writeup/data/riemann/lzeros'
def theta(t, q, a): return np.imag(loggamma(0.25 + a / 2 + 0.5j * t)) + t / 2 * math.log(q / math.pi)
def check(g, q, a, lo, hi, label):
    g = np.sort(g); gb = g[(g >= lo) & (g <= hi)]
    u = gb / (2 * math.pi) * np.log(q * gb / (2 * math.pi * math.e)); s = np.diff(u); o = np.argsort(s)[::-1][:4]
    t = np.linspace(lo, hi, 9); drift = np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi
    print(label, 'zeros in band %d; largest unfolded gaps %s at heights %s; N - theta/pi at %s: %s' % (len(gb),
          ' '.join('%.2f' % s[i] for i in o), ' '.join('%.1f' % gb[i] for i in o), list(t.astype(int)), ' '.join('%+.2f' % d for d in drift)))
Z = dict(np.load(LZ + '/r134_C.npz'))
for D in (-195, 195):
    if 'D%d' % D in Z: check(Z['D%d' % D], 195, 0 if D > 0 else 1, 600, 1000, 'C  D=%d' % D)
for D in (145, -145, 247, -247):
    if 'D%d' % D in Z: check(Z['D%d' % D], abs(D), 0 if D > 0 else 1, 600, 1000, 'C  D=%d (control)' % D)
Z = dict(np.load(LZ + '/r136_Dn.npz'))
for D in (939, -939):
    if 'D%d' % D in Z: check(Z['D%d' % D], 939, 0 if D > 0 else 1, 100, 600, 'Dn D=%d' % D)
Z = {}
for i in range(4):
    z = np.load(LZ + '/r121_part%d.npz' % i); Z.update({k: z[k] for k in z.files})
for k in Z:
    if k.startswith('N173_'):
        for side in (1, -1):
            g = side * Z[k][side * Z[k] > 0]
            u = np.sort(g); gb = u[(u >= 100) & (u <= 500)]; uu = gb / (2 * math.pi) * np.log(173 * gb / (2 * math.pi * math.e)); s = np.diff(uu); o = np.argsort(s)[::-1][:4]
            print('A  %s side %+d: zeros in band %d; largest unfolded gaps %s at %s' % (k, side, len(gb), ' '.join('%.2f' % s[i] for i in o), ' '.join('%.1f' % gb[i] for i in o)))
