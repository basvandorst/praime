"""Expected one-point readings from the primes alone: the frozen resonance sums R put in place of the zero sums G
(explicit formula; no zero used), and under the alternatives (no kappa: -T/2pi; no root number: -kappa T/2pi)."""
import sys, os, math, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import swap2_test as S
F = S.load_frozen(HERE); O = F['onepoint.npz']; out = {}
for des in ('OO', 'HI', 'C'):
    funcs = S.design_functions(des, F); key = des
    idx = {str(l): k for k, l in enumerate(O[key + '_labels'])}
    G = np.zeros((len(funcs), 2), complex)
    for k, f in enumerate(funcs):
        if f['label'] in idx:
            i = idx[f['label']]; G[k] = (O[key + '_res_sharp'][i], O[key + '_res_hann'][i])
    for hyp in ('primes', 'law'):
        GG = G.copy()
        if hyp == 'law':
            for k, f in enumerate(funcs):
                if f['label'] in idx:
                    i = idx[f['label']]; GG[k] = (O[key + '_law_sharp'][i], O[key + '_law_hann'][i])
        r = S.readings_onepoint(des, F, GG, np.zeros(len(funcs), bool), funcs)
        print(des, hyp)
        for k, v in r.items(): print(S.fmt(k, v))
        out['%s_%s' % (des, hyp)] = {k: S.jsonable(v) for k, v in r.items()}
json.dump(out, open(os.path.join(HERE, 'onepoint_expected.json'), 'w'), indent=1)
