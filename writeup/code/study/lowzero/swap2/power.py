"""Expected value of every registered pair reading under each hypothesis (noise-free, from the frozen predictions) and its
expected standard error (noise model fitted on the USED set Fn: per-function additive and multiplicative parts plus a
pair part, Re and Im separately, scaled by sqrt(L_i L_j / L_Fn^2 * 500/T)). No reserved zero is used.
usage: power.py DESIGN NREP   -> power_DESIGN.json"""
import sys, os, math, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import swap2_test as S

des = sys.argv[1]; nrep = int(sys.argv[2]) if len(sys.argv) > 2 else 20
F = S.load_frozen(HERE); funcs = S.design_functions(des, F); n = len(funcs)
NM = json.load(open(os.path.join(HERE, 'noise_model.json')))
W1, W2 = S.DESIGN[des]['win']; T = W2 - W1; tm = (W1 + W2) / 2
Lf = np.array([math.log(f['q'] * tm / (2 * math.pi)) for f in funcs]); LFN = math.log(1150 * 350 / (2 * math.pi))
excl = np.zeros(n, bool)
t0 = time.time()


def pairs(sig):
    """(I, J, {hyp: values}) for the design"""
    if des in ('OO', 'HI', 'FN'):
        P = F['pred_FN.npz' if des == 'FN' else 'pred_OO.npz']; I, J = P['I'], P['J']
        m = (I < 181) & (J < 181) if des == 'HI' else np.ones(len(I), bool)
        p = P['pred_%.1f' % sig][m]; p1 = P['pred1_%.1f' % sig][m]
        return I[m], J[m], {'law': p, 'half': 0.5 * p, 'none': 0 * p, 'const': p1}
    P = F['pred_C.npz']
    I = np.concatenate([P['I_cc'], P['I_cr']]); J = np.concatenate([P['J_cc'], P['J_cr']])
    g = lambda k: np.concatenate([P['%s_cc_%.1f' % (k, sig)], P['%s_cr_%.1f' % (k, sig)]])
    p = g('pred')
    return I, J, {'law': p, 'half': 0.5 * p, 'none': 0 * p, 'const': g('pred1'), 'noroot': g('prednoroot'), 'conjJ': g('predconjJ')}


def make_X(vals_by_sig, rng=None):
    X = {}
    for sig in S.SIGMAS:
        I, J, val = vals_by_sig[sig]
        M = np.zeros((n, n), complex)
        v = val.astype(complex).copy()
        if rng is not None:
            nm = NM['%.1f' % sig]; f = np.sqrt(Lf[I] * Lf[J] / LFN ** 2 * 500 / T)
            comp = []
            for part in ('re', 'im'):
                c = nm[part]; u = rng.normal(0, c['sd_u2'], n); w = rng.normal(0, c['sd_v'], n)
                comp.append(u[I] + u[J] + rng.normal(0, c['sd_e2'], len(I)))
                comp.append(w[I] + w[J])
            mult = comp[1] + 1j * comp[3]
            v = v + f * (comp[0] + 1j * comp[2] + mult * vals_by_sig['law_%.1f' % sig])
        M[I, J] = v; M[J, I] = np.conj(v)
        X[sig] = M
    return X


hyps = None; base = {}
for sig in S.SIGMAS:
    I, J, hv = pairs(sig); base[sig] = (I, J, hv); hyps = list(hv)
out = {'design': des, 'expected': {}, 'se': {}}
for h in hyps:
    vs = {sig: (base[sig][0], base[sig][1], base[sig][2][h]) for sig in S.SIGMAS}
    for sig in S.SIGMAS: vs['law_%.1f' % sig] = base[sig][2]['law']
    r = S.readings_pairs(des, F, make_X(vs), excl, funcs)
    out['expected'][h] = {k: [[float(complex(c).real), float(complex(c).imag)] for c in np.atleast_1d(v['coef'])] for k, v in r.items()}
    print('%s: expected under %s done (%.0fs)' % (des, h, time.time() - t0), flush=True)
rng = np.random.default_rng(42)
sims = {}; jks = {}
vs = {sig: (base[sig][0], base[sig][1], base[sig][2]['law']) for sig in S.SIGMAS}
for sig in S.SIGMAS: vs['law_%.1f' % sig] = base[sig][2]['law']
for rep in range(nrep):
    r = S.readings_pairs(des, F, make_X(vs, rng), excl, funcs)
    for k, v in r.items():
        sims.setdefault(k, []).append(np.atleast_1d(v['coef'])); jks.setdefault(k, []).append(np.concatenate([v['se_re'], v['se_im']]))
    print('  rep %d (%.0fs)' % (rep, time.time() - t0), flush=True)
for k in sims:
    a = np.array(sims[k]); j = np.array(jks[k])
    out['se'][k] = {'sim_sd_re': a.real.std(0).tolist(), 'sim_sd_im': a.imag.std(0).tolist(), 'mean_re': a.real.mean(0).tolist(),
                    'mean_im': a.imag.mean(0).tolist(), 'jack_mean': j.mean(0).tolist()}
json.dump(out, open(os.path.join(HERE, 'power_%s.json' % des), 'w'), indent=1)
# compact table
print('reading | expected under ' + ' | '.join(hyps) + ' || se (sim sd of the last coefficient, mean jackknife se)')
for k in out['expected'][hyps[0]]:
    ex = ' | '.join(('%.3f' % out['expected'][h][k][-1][0]) + (('%+.3fi' % out['expected'][h][k][-1][1]) if des == 'C' else '') for h in hyps)
    se = out['se'][k]
    print('%-66s %s || %.4f, %.4f' % (k, ex, se['sim_sd_re'][-1], se['jack_mean'][len(se['sim_sd_re']) - 1]))
