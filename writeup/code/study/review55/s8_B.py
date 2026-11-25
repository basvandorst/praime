"""B_good = prod_{good p <= PB} F_p(1/p) L_p(1, pi x pi~)  (absolutely convergent; bad p set to 1) for all degree >= 2 objects;
alternative constant 1/L(1, Ad) ~ c / B.  Then chi^2 of the measured c_hat against the law, c = 1, c = 1/L(1,Ad) (no B), c = 0,
for degree >= 2 and >= 3, with the registered sd and with the sd rescaled by the measured Var/M per family."""
import json, math, os, sys
import numpy as np, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s6_const import F_of, satake_sym
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'; D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json'))
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
def B_of(alpha, ps):
    X = 1 / ps.astype(float); Fp = F_of(alpha, X).real
    LRS = np.ones(len(ps))
    d = alpha.shape[1]
    for i in range(d):
        for j in range(d):
            LRS = LRS / (1 - (alpha[:, i] * np.conj(alpha[:, j])).real * X)  # real for self-dual sets (pairs conj)
    return float(np.prod(Fp * LRS))
out = {}
cache = D + '/s8_B.json'
if os.path.exists(cache): out = json.load(open(cache))
for key, e in P['L'].items():
    if e['degree'] < 2 or key in out: continue
    sp = e['spec']
    if sp[0] == 'mf':
        PB = 2000; k, N, idx = sp[1], sp[2], sp[3]
        mf = pari.mfinit([N, k], 0); co = pari.mfcoefs(pari.mfeigenbasis(mf)[idx], PB)
        ps = np.array([int(x) for x in pari('primes([2,%d])' % PB)]); ps = ps[N % ps != 0]
        lam = np.array([float(co[int(p)]) for p in ps]) / ps.astype(float) ** ((k - 1) / 2)
        ph = np.arccos(np.clip(lam / 2, -1, 1)); alpha = np.stack([np.exp(1j * ph), np.exp(-1j * ph)], 1)
    elif sp[0] == 'sym':
        PB = 2000; E = pari.ellinit(sp[1]); N = int(pari.ellglobalred(E)[0]); m = sp[2]
        an = np.array([int(x) for x in pari.ellan(E, PB)], float)
        ps = np.array([int(x) for x in pari('primes([2,%d])' % PB)]); ps = ps[N % ps != 0]
        alpha = satake_sym(an[ps - 1], ps.astype(float), m)
    else:
        PB = 1000; N = int(e['Q']); Fp = 4 * pari(sp[1]) + pari(sp[2]) ** 2
        ps = np.array([int(x) for x in pari('primes([3,%d])' % PB)]); ps = ps[N % ps != 0]
        alpha = np.array([np.roots([float(x) for x in pari.Vec(pari.hyperellcharpoly(Fp * pari.Mod(1, int(p))))]) / math.sqrt(p) for p in ps])
    out[key] = dict(B=B_of(alpha, ps), PB=PB)
    json.dump(out, open(cache, 'w'))
R = json.load(open(D + '/s1_rows.json'))
var_along = {'T2': 0.192, 'T3': 0.330, 'T4a': 0.479, 'T4b': 0.152}     # measured (s3b_var.py)
def chi(keys, scale):
    Y = []; W = []; C = []; A = []
    for k in keys:
        r = R[k]; eps = complex(*r['eps']); unit = r['T2'] / (2 * math.pi)
        y = (complex(*r['G']) * (-eps).conjugate()).real / unit
        s = r['sd_along'] / unit * (math.sqrt(var_along[r['cls']] / 0.077) if scale else 1.0)
        Y.append(y); W.append(1 / s ** 2); C.append(r['c']); A.append(r['c'] / out[k]['B'])
    Y, W, C, A = map(np.array, (Y, W, C, A))
    f = lambda pred: float(np.sum(W * (Y - pred) ** 2))
    # free amplitude fits: y = b * model
    fit = lambda m: (float(np.sum(W * Y * m) / np.sum(W * m * m)), f(np.sum(W * Y * m) / np.sum(W * m * m) * m))
    return dict(n=len(keys), law=f(C), c1=f(np.ones(len(C))), noB=f(A), zero=f(np.zeros(len(C))), fit_law=fit(C), fit_noB=fit(A), fit_c1=fit(np.ones(len(C))),
                Brange=(float(min(out[k]['B'] for k in keys)), float(max(out[k]['B'] for k in keys))), corr_c_A=float(np.corrcoef(C, A)[0, 1]))
reg = [k for k in R if R[k]['status'] == 'COMPLETE']
for name, ks in (('degree >= 2', [k for k in reg if R[k]['deg'] >= 2]), ('degree >= 3', [k for k in reg if R[k]['deg'] >= 3]),
                 ('degree 2', [k for k in reg if R[k]['deg'] == 2]), ('degree 3', [k for k in reg if R[k]['deg'] == 3]), ('degree 4', [k for k in reg if R[k]['deg'] == 4])):
    for sc in (False, True):
        d = chi(ks, sc)
        print('%-12s n=%2d %-14s chi2: law %7.1f | c=1 %7.1f | 1/L(1,Ad) no B %7.1f | c=0 %7.1f || free-amplitude fits: law b=%.3f chi2 %.1f; noB b=%.3f chi2 %.1f; c=1 b=%.3f chi2 %.1f | B in [%.3f, %.3f], corr(c, c/B) %.3f' % (
            name, d['n'], 'measured sd' if sc else 'registered sd', d['law'], d['c1'], d['noB'], d['zero'], d['fit_law'][0], d['fit_law'][1], d['fit_noB'][0], d['fit_noB'][1], d['fit_c1'][0], d['fit_c1'][1], *d['Brange'], d['corr_c_A']))
