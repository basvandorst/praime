"""Review 31: window averages over ln X of R(X, d) from the zeros: diagonal (=D, closed form) + window mean of
the off-diagonal sum (2/d) Re sum_{j!=k} c_j conj(c_k) X^{i(g_j-g_k)}/(2+i(g_j-g_k)), averaged exactly over u = ln X in [a, b]."""
import sys, math, numpy as np
nz = int(sys.argv[1])
Z = np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt') if l.strip()])[:nz]
G = np.concatenate([Z, -Z]); RHO = 0.5 + 1j * G
def closed(d):
    return math.log(1/d) - 0.5772156649015329 - math.log(2*math.pi) + 2*math.log1p(d) + ((2+d)/d)*math.log1p(d/2)
def offwin(d, a, b):
    c = ((1 + d) ** RHO - 1) / RHO; L = b - a; tot = 0.0
    ea = np.exp(1j * G * a); eb = np.exp(1j * G * b)
    for i0 in range(0, len(G), 400):
        Dg = G[i0:i0+400, None] - G[None, :]
        M = c[i0:i0+400, None] * np.conj(c)[None, :] / (2 + 1j * Dg)
        num = eb[i0:i0+400, None] * np.conj(eb)[None, :] - ea[i0:i0+400, None] * np.conj(ea)[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            W = num / (1j * Dg * L)
        W[Dg == 0] = 0.0
        tot += float(np.real(np.sum(M * W)))
    return 2 * tot / d
l7, l10 = math.log(1e7), math.log(1e10)
for d in (0.1, 0.01):
    print('delta', d, 'D', round(closed(d), 4), flush=True)
    for (a, b, lab) in [(l7, l10, '[1e7,1e10]')] + [(l7 + k*(l10-l7), l7 + (k+1)*(l10-l7), 'win%d' % k) for k in range(1, 12)] + [(l7, math.log(1e40), '[1e7,1e40]'), (l7, math.log(1e100), '[1e7,1e100]')]:
        print('  %-12s off-diag window mean %+.4f' % (lab, offwin(d, a, b)), flush=True)
