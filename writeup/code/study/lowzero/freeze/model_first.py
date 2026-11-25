"""model_first.py: per-prime contrasts (p missing minus p present, the other small primes present) of the interval
observables in the random models at one size; levels of the all-present field; and the closed-form alternatives.
usage: model_first.py N B seed"""
import sys, os, json, time
import numpy as np
from scipy.special import hyp2f1
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modellib as ml

PR = (2, 3, 5, 7)
KEYS = ['M', 'lnZ_0.5', 'lnZ_1', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'mean_ln']


def explicit_primes(k0):
    return [p for p in ml.primes_upto(200) if np.log(p) < k0 - 0.5]


def run(model, N, B, rng):
    L = float(N)
    if model.startswith('hyb'):
        k0 = int(model[3:])
        W = ml.cue_field(rng, B, N, k0=k0)
        ex = explicit_primes(k0)
    elif model == 'cue':
        W = ml.cue_field(rng, B, N, k0=1); ex = list(PR)
    elif model == 'ep':
        W = ml.ep_field(rng, B, L, skip=PR); ex = list(PR)
    Y = {}
    for p in ex:
        Y[p] = ml.probe(rng, B, p)
        W = W + Y[p]
    s_all = ml.stats(W)
    res = {'level': {k: [float(s_all[k].mean()), float(s_all[k].std())] for k in KEYS}}
    for p in [q for q in PR if q in Y]:
        s_m = ml.stats(W - Y[p])
        res[p] = {}
        for k in KEYS:
            d = s_m[k] - s_all[k]
            res[p][k] = [float(d.mean()), float(d.std() / np.sqrt(B))]
        res[p]['varM'] = float(s_m['M'].var() - s_all['M'].var())
    return res


def closed_forms():
    out = {}
    phis = (np.arange(2000) + 0.5) * 2 * np.pi / 2000
    h = (np.arange(4000) + 0.5) * ml.THETA / 4000
    for p in PR:
        Yp = -np.log(np.abs(1 - p ** -0.5 * np.exp(1j * (phis[:, None] + h[None, :] * np.log(p)))))
        o = {}
        for b in ml.BETAS:
            o['naive_lnZ_%g' % b] = float(-np.log(hyp2f1(b, b, 1, 1.0 / p)))
            o['spread_lnZ_%g' % b] = float(-np.mean(np.log(np.mean(np.exp(2 * b * Yp), axis=1))))
            Kp2 = -np.log(1 - 1.0 / p)
            o['brw_lnZ_%g' % b] = float(-np.log(hyp2f1(b, b, 1, 1.0 / p))) if b <= 1 else float(-b * Kp2)
        o['brw_M'] = float(0.5 * np.log(1 - 1.0 / p))
        o['spread_M'] = float(-np.mean(Yp.max(axis=1)))
        o['var_Yp'] = float(Yp.var())
        o['var_intervalmean_Yp'] = float(Yp.mean(axis=1).var())
        out[p] = o
    return out


def main():
    N = int(sys.argv[1]); B = int(sys.argv[2]); seed = int(sys.argv[3])
    rng = np.random.default_rng(seed)
    cf = closed_forms()
    print('closed forms (missing minus present):')
    for p in PR:
        print(p, ' '.join('%s %.4f' % (k, v) for k, v in cf[p].items()))
    allres = {'closed': cf}
    for model in ['hyb2', 'hyb3', 'hyb4', 'cue', 'ep']:
        t = time.time()
        r = run(model, N, B, rng)
        allres[model] = r
        print('\n%s N=%d B=%d (%.1f s)  levels: %s' % (model, N, B, time.time() - t,
              ' '.join('%s %.3f/%.3f' % (k, *r['level'][k]) for k in KEYS)))
        for p in [q for q in PR if q in r]:
            print('  p=%d ' % p + ' '.join('%s %+.4f(%.4f)' % (k, *r[p][k]) for k in KEYS) + ' varM %+.4f' % r[p]['varM'])
    json.dump(allres, open('model_first_N%d.json' % N, 'w'), indent=1)


if __name__ == '__main__':
    main()
