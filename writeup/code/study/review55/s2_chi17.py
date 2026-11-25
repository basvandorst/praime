"""The two excluded characters mod 17: recompute zeros with a finer subdivision, compute arg L(1/2, chi) by continuous
variation along the real axis from sigma = 40, and give the predicted offset 1/2 - A0/pi.  Then the law's sums."""
import json, math, sys, os, time
import numpy as np, cypari2
from scipy.special import loggamma
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'
P = json.load(open(B + '/predictions_blind.json'))
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6); pari.set_real_precision(38)
def theta(t, Q, mu):
    t = np.asarray(t, float); s = t / 2 * math.log(Q / math.pi ** len(mu))
    for m in mu: s = s + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return s
keys = sys.argv[1:] or ['T1_17_5', 'T1_17_7', 'T1_23_7']
for key in keys:
    e = P['L'][key]; q, a = e['spec'][1], e['spec'][2]
    L = pari('lfuncreate([G=znstar(%d,1), znconreychar(G,%d)])' % (q, a))
    eps = complex(pari.lfunrootres(L)[2])
    # continuous arg of L(sigma) from sigma = 40 down to 1/2
    sig = np.linspace(40, 0.5, 4001); args = []
    for s in sig:
        v = complex(pari.lfun(L, s)); args.append(math.atan2(v.imag, v.real))
    A = np.unwrap(np.array(args)); A0 = A[-1]
    v12 = complex(pari.lfun(L, 0.5))
    c_true = 0.5 - A0 / math.pi
    t0 = time.time()
    Li = pari.lfuninit(L, [1, 1003])
    gf = np.sort(np.array([float(x) for x in pari.lfunzeros(Li, [0.0001, 1001.0], 32)]))
    secs = time.time() - t0
    z = np.load(B + '/zeros_blind/' + key + '.npz'); g = z['g']
    # match
    n = min(len(g), len(gf))
    same = len(g) == len(gf) and np.max(np.abs(g - gf)) < 1e-6
    Q, mu, T2 = e['Q'], e['mu'], e['T2']
    r = np.arange(1, len(gf) + 1) - theta(gf, Q, mu) / math.pi
    gg = gf[(gf > 0) & (gf <= T2)]
    ph = np.exp(2j * theta(gg, Q, mu)); w = np.sin(math.pi * gg / T2) ** 2
    G = complex(ph.sum()); Gh = complex((w * ph).sum()); pred = complex(*e['pred'])
    print('%s: eps %.4f%+.4fi  arg eps %.4f  L(1/2) = %.6f%+.6fi  A0 (continuous) = %.4f  -> c_exp = %.4f (registered rule %.4f)' % (
        key, eps.real, eps.imag, math.atan2(eps.imag, eps.real), v12.real, v12.imag, A0, c_true, json.loads(str(z['info']))['c_exp']))
    print('   fine search divz 32 on (1e-4, 1001]: %d zeros (file %d), identical to 1e-6: %s, first zero %.6f, %.0f s; r_mean %.4f' % (
        len(gf), len(g), same, gf[0], secs, r.mean()))
    print('   law: G = %+.2f%+.2fi  pred %+.2f%+.2fi  slope Re(G/pred) = %.4f ± %.4f  Hann slope %.4f   im %.4f' % (
        G.real, G.imag, pred.real, pred.imag, (G / pred).real, e['sd']['along'] / abs(pred), (Gh / (pred / 2)).real, (G / pred).imag))
    json.dump(dict(G=[G.real, G.imag], Gh=[Gh.real, Gh.imag], n=len(gg), A0=A0, c_exp=c_true), open(os.path.dirname(os.path.abspath(__file__)) + '/s2_%s.json' % key, 'w'))
