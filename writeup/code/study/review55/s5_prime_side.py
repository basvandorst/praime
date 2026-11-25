"""Prime side of the explicit formula (fact 2056 (i)) for the Artin L-function of a cubic field, by stationary phase:
   sum_gamma H(gamma) e^{2i theta(gamma)} = -(1/2pi) sum_n Lambda(n) n^{-1/2} int H(t) e^{i(2theta(t) - t ln n)} dt + E0,
   int ~ H(t_n) sqrt(2pi/phi'') e^{i(phi(t_n) + pi/4)},  2theta'(t_n) = ln n.   Lambda_rho(p^k) = ln p (sum_{f_i | k} f_i - 1).
Compared field by field with the zero side (rho zeros = zeta_K zeros minus zeta zeros) for several windows inside (0, 300]."""
import json, math, os, sys
import numpy as np, cypari2
from scipy.special import loggamma
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'; D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json'))
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
zz = np.load(D + '/zeta_zeros_310.npy')
def theta(t, Q, mu):
    t = np.asarray(t, float); s = t / 2 * math.log(Q / math.pi ** len(mu))
    for m in mu: s = s + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return s
def th1(t, Q, mu, h=1e-4): return (theta(t + h, Q, mu) - theta(t - h, Q, mu)) / (2 * h)
def th2(t, Q, mu, h=1e-3): return (theta(t + h, Q, mu) - 2 * theta(t, Q, mu) + theta(t - h, Q, mu)) / h ** 2
T2 = 300.0
WIN = {
    'hann':  lambda t: np.where((t > 0) & (t < T2), np.sin(np.pi * t / T2) ** 2, 0.0),
    'low':   lambda t: np.where((t > 0) & (t < 100), np.sin(np.pi * t / 100) ** 2, 0.0),
    'mid':   lambda t: np.where((t > 100) & (t < 200), np.sin(np.pi * (t - 100) / 100) ** 2, 0.0),
    'high':  lambda t: np.where((t > 200) & (t < 300), np.sin(np.pi * (t - 200) / 100) ** 2, 0.0),
    'tukey': lambda t: np.where((t > 0) & (t < T2), np.where(t < 40, np.sin(np.pi * t / 80) ** 2, np.where(t > T2 - 40, np.sin(np.pi * (T2 - t) / 80) ** 2, 1.0)), 0.0),
}
INT = {'hann': T2 / 2, 'low': 50.0, 'mid': 50.0, 'high': 50.0, 'tukey': T2 - 40.0}
def lam_rho(polstr, d, nmax):
    """dict-free arrays: n, Lambda(n), is_prime for prime powers n <= nmax."""
    pol = pari(polstr); disc = int(pari.poldisc(pol)); nf = pari.nfinit(pol)
    ps = np.array([int(p) for p in pari('primes([2,%d])' % nmax)], dtype=np.int64)
    fs = []
    for p in ps:
        p = int(p)
        if disc % p == 0:
            dec = pari.idealprimedec(nf, p); fs.append(tuple(int(pr[3]) for pr in dec))
        else:
            r = len(pari.polrootsmod(pol, p))
            fs.append((1, 1, 1) if r == 3 else ((1, 2) if r == 1 else (3,)))
    N = []; L = []; K = []
    for p, f in zip(ps, fs):
        p = int(p); lp = math.log(p); pk = p; k = 1
        while pk <= nmax:
            N.append(pk); L.append(lp * (sum(fi for fi in f if k % fi == 0) - 1)); K.append(k)
            pk *= p; k += 1
    return np.array(N, float), np.array(L), np.array(K), fs, ps
def prime_side(Q, mu, N, L, K):
    # stationary points by Newton, start 2 pi sqrt(n/Q) (degree 2)
    t = 2 * np.pi * np.sqrt(N / Q)
    for _ in range(30):
        f = 2 * th1(t, Q, mu) - np.log(N); fp = 2 * th2(t, Q, mu)
        t = np.maximum(t - f / fp, 1e-3)
    phi2 = 2 * th2(t, Q, mu)
    base = -(1 / (2 * np.pi)) * L * N ** -0.5 * np.sqrt(2 * np.pi / phi2) * np.exp(1j * (2 * theta(t, Q, mu) - t * np.log(N) + np.pi / 4))
    res = {}
    for w, H in WIN.items():
        h = H(t)
        res[w] = dict(all=complex((base * h).sum()), primes=complex((base * h)[K == 1].sum()), squares=complex((base * h)[K == 2].sum()),
                      higher=complex((base * h)[K >= 3].sum()))
    return res, float(np.abs(2 * th1(t, Q, mu) - np.log(N)).max())
keys = [k for k, e in P['products'].items() if e['cls'] == 'T5' and e['rho']['c'] != 0]
if len(sys.argv) > 1: keys = [k for k in keys if k in sys.argv[1:]]
out = {}
for key in keys:
    e = P['products'][key]; r = e['rho']; Q, mu, c = r['Q'], r['mu'], r['c']
    nmax = int(math.exp(2 * th1(np.array([T2 + 2.0]), Q, mu)[0])) + 1
    N, L, K, fs, ps = lam_rho(e['spec'][1], e['d'], nmax)
    ps_res, newton_err = prime_side(Q, mu, N, L, K)
    g = np.load(B + '/zeros_blind/' + key + '.npz')['g']
    isz = np.array([np.min(np.abs(zz - x)) < 1e-6 for x in g]); gr = g[~isz]; gr = gr[(gr > 0) & (gr <= T2)]
    ph = np.exp(2j * theta(gr, Q, mu))
    row = dict(d=e['d'], c=c, nmax=nmax, newton=newton_err)
    for w, H in WIN.items():
        Z = complex((H(gr) * ph).sum()); law = -c * INT[w] / (2 * np.pi)
        row[w] = dict(zero=[Z.real, Z.imag], prime=[ps_res[w]['all'].real, ps_res[w]['all'].imag], law=law,
                      primes_only=ps_res[w]['primes'].real, squares=ps_res[w]['squares'].real, higher=ps_res[w]['higher'].real)
    # split-type counts for p <= 1000 (check of the splitting law 1/6, 1/2, 1/3)
    sm = [f for p, f in zip(ps, fs) if p < 100000 and e['d'] % int(p) != 0]
    row['split_frac'] = [sum(1 for f in sm if f == (1, 1, 1)) / len(sm), sum(1 for f in sm if f == (1, 2)) / len(sm), sum(1 for f in sm if f == (3,)) / len(sm)]
    out[key] = row
    print('%-9s d=%4d c=%.3f nmax=%7d  ' % (key, e['d'], c, nmax) + '  '.join(
        '%s Z %+6.1f%+5.1fi P %+6.1f%+5.1fi law %+6.1f' % (w, row[w]['zero'][0], row[w]['zero'][1], row[w]['prime'][0], row[w]['prime'][1], row[w]['law'])
        for w in ('hann', 'low', 'high')), flush=True)
json.dump(out, open(D + '/s5_rows%s.json' % ('' if len(sys.argv) == 1 else '_' + '_'.join(sys.argv[1:])), 'w'), indent=0)
ks = list(out)
print('\npooled over %d fields (slope on the law = sum X*law / sum law^2):' % len(ks))
for w in WIN:
    den = sum(out[k][w]['law'] ** 2 for k in ks)
    sz = sum(out[k][w]['zero'][0] * out[k][w]['law'] for k in ks) / den
    sp = sum(out[k][w]['prime'][0] * out[k][w]['law'] for k in ks) / den
    s1 = sum(out[k][w]['primes_only'] * out[k][w]['law'] for k in ks) / den
    s2 = sum(out[k][w]['squares'] * out[k][w]['law'] for k in ks) / den
    s3 = sum(out[k][w]['higher'] * out[k][w]['law'] for k in ks) / den
    rms = math.sqrt(np.mean([(out[k][w]['zero'][0] - out[k][w]['prime'][0]) ** 2 + (out[k][w]['zero'][1] - out[k][w]['prime'][1]) ** 2 for k in ks]))
    print('%-6s zero side %.4f  prime side %.4f  (primes %.4f, squares %+.4f, higher powers %+.4f)  rms |Z - P| per field %.2f' % (w, sz, sp, s1, s2, s3, rms))
