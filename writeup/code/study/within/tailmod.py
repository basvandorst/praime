"""Long-range tail of the within part from insertion responses (linear response, HL pair correlations of the primes > y),
for any sieve level y (chunked; no kernel).  See chainlaw.py for the formula."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math
from hlcore import neighbour_class_laws, S_tail_table
from chainlaw import proj_within, increment_vectors
from scipy.stats import poisson
def responses_all(y, lam, qs, T=700):
    Q, U, PRf, PLf, rho0 = neighbour_class_laws(y, lam, qs, T=T, chunk=2048)
    out = {}
    for q in qs:
        red = [a for a in range(q) if math.gcd(a, q) == 1]; C = len(red)
        PR = PRf[q][:, red]; PL = PLf[q][:, red]; lut = -np.ones(q, int); lut[red] = np.arange(C); cu = lut[U % q]
        Pm = np.array([PR[cu == a].mean(0) for a in range(C)]); Pm /= Pm.sum(1, keepdims=True); pi = np.full(C, 1 / C)
        PiW = proj_within(Pm, pi); w, Ub = np.linalg.eigh(PiW); Ub = Ub[:, w > 0.5]
        Fv = (increment_vectors(Pm, pi) @ Ub).reshape(C, C, -1)
        r = np.einsum('ua,uad->ud', PL, Fv[:, cu].transpose(1, 0, 2)) + np.einsum('ub,ubd->ud', PR, Fv[cu]) - np.einsum('ua,ub,abd->ud', PL, PR, Fv)
        out[q] = dict(r=r, C=C, Pm=Pm)
    return Q, U, out
def autocorr(Q, U, r):
    A = np.zeros(Q)
    for d in range(r.shape[1]):
        full = np.zeros(Q); full[U] = r[:, d]; F = np.fft.rfft(full); A += np.fft.irfft(F * np.conj(F), n=Q)
    return A / len(U)
def tail_values(Q, U, A, C, lam, y, Hs, lagbins=None, weight='w2', Stab=None):
    pref = (Q / len(U)) / ((C - 1) ** 2 * lam)
    Hmax = max(Hs); h = np.arange(1, Hmax + 1)
    if Stab is None: Stab = S_tail_table(Hmax, y)
    mu = h / lam
    w = 1 - np.exp(-mu) if weight == 'w2' else 1 - np.exp(-mu) * (1 + mu)
    base = A[h % Q] * (Stab[1:Hmax + 1] - 1)
    res = {'V_per_dir': float(A[0]) / (C - 1) ** 2,
           'W_tail': {H: float(2 * pref * np.sum((1 - h[:H] / H) * w[:H] * base[:H])) for H in Hs}}
    if lagbins:
        res['bins'] = {}
        for lo, hi in lagbins:
            wg = poisson.cdf(hi - 2, mu) - poisson.cdf(lo - 2, mu)
            res['bins'][(lo, hi)] = float(pref * np.sum(wg * base))
    return res
if __name__ == '__main__':
    import json, time
    lam = math.log(1.5e11 + 5e7)
    BINS = [(3, 6), (6, 11), (11, 31), (31, 101), (101, 301), (301, 1001), (1001, 3001)]
    tl = json.load(open('data/val_tail.json'))
    Hs = [10 ** 4, 10 ** 5, 10 ** 6]
    for y in (11, 13, 17):
        t0 = time.time(); Q, U, R = responses_all(y, lam, (3, 5, 7)); Stab = S_tail_table(10 ** 6, y)
        for q in (3, 5, 7):
            A = autocorr(Q, U, R[q]['r'])
            tv = tail_values(Q, U, A, R[q]['C'], lam, y, Hs, BINS, Stab=Stab)
            tv3 = tail_values(Q, U, A, R[q]['C'], lam, y, Hs, None, weight='w3', Stab=Stab)
            print('y=%2d q=%d V/dir %.4f  W_tail(1e4,1e5,1e6) w2: %+.4f %+.4f %+.4f   w3: %+.4f %+.4f %+.4f' % (y, q, tv['V_per_dir'], *[tv['W_tail'][H] for H in Hs], *[tv3['W_tail'][H] for H in Hs]))
            mb = dict(zip([tuple(b) for b in tl['bins']], tl['%d_1' % q]['binsum']))
            print('      bins model/measured: ' + '  '.join('[%d,%d) %+.4f/%+.4f' % (lo, hi, tv['bins'][(lo, hi)], mb[(lo, hi)][0]) for lo, hi in BINS))
        print('   (%.1fs)' % (time.time() - t0), flush=True)
