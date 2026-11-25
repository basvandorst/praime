# Alternative amplitude models on all data: 20 primitive L (this test) + 138 curves (round 148, band (20,400]).
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, math, numpy as np
from scipy.special import loggamma
D = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/'
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
r2 = json.load(open(OUT + 'r2_noise.json'))
r4 = json.load(open(OUT + 'r4_p2.json'))
rows = []
for key, v in r2.items():
    P = pred[key]; cj = json.load(open(D + 'constants_%s.json' % key))
    A1 = cj.get('A1'); 
    if v['fam'].startswith('GL1'): A1 = 1.0
    se = np.sqrt(np.mean([v['jk'][k][0] ** 2 for k in ('4', '6', '8', '12', '20', '30')]))
    var = (se * v['absp']) ** 2
    G = complex(*v['G'])
    P2 = complex(*r4[key]['P2']) if key in r4 else 0j
    rows.append(dict(key=key, fam=v['fam'], G=G, Gc=G - P2, W=P['T2'] / (2 * math.pi), eps=complex(*P['eps']), c=P['c'], A1=A1, var=var, M=v['M']))
    print('%-5s c=%.4f A1=%s  se_K-avg=%.3f' % (key, P['c'], A1, se))
# round 148
F = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r148/fresh_onepoint.npz')
Pf = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/ecross/pred_fresh.npz')
Gl = F['Gl'][:, 0, 0]
W = 380 / (2 * math.pi)
def thE(N, t): return t / 2 * math.log(N) - t * math.log(2 * math.pi) + np.imag(loggamma(1 + 1j * t))
lawE = Pf['lineamp'][:, 0] * W
assert np.allclose(Pf['lineamp'][:, 0], -Pf['w'] * Pf['c'])
ME = np.array([(thE(N, 400.0) - thE(N, 20.0)) / math.pi for N in Pf['N']])
resE = Gl.real - lawE
kapE = np.mean(resE ** 2 / ME)
print('round 148: slope sum(Re G law)/sum law^2 = %.4f ; kappa (along) from residuals = %.3f' % ((Gl.real * lawE).sum() / (lawE ** 2).sum(), kapE))
for i in range(len(Gl)):
    rows.append(dict(key=str(Pf['ids'][i]), fam='EC r148', G=complex(Gl[i]), Gc=complex(Gl[i]), W=W, eps=complex(Pf['w'][i]), c=float(Pf['c'][i]),
                     A1=float(Pf['LS1'][i]), var=kapE * ME[i], M=ME[i]))

def chi2(model, rowsel, use_corr=False):
    s = 0.0; n = 0
    for r in rowsel:
        A = model(r)
        if A is None: continue
        G = r['Gc'] if use_corr else r['G']
        # complex residual: along+across, each with variance var (use var for both components)
        res = G - A * r['W']
        s += (res.real ** 2 + res.imag ** 2) / r['var']; n += 2
    return s, n

models = {
 'law  -eps c': lambda r: -r['eps'] * r['c'],
 'c = 1 for all: -eps': lambda r: -r['eps'],
 'no B: -eps/L(1,Ad)': lambda r: (-r['eps'] / r['A1']) if r['A1'] else None,
 'no root number: -c': lambda r: -r['c'],
 'conj eps: -conj(eps) c': lambda r: -np.conj(r['eps']) * r['c'],
 'H0: 0': lambda r: 0.0,
}
for sel_name, sel in (('20 primitive', [r for r in rows if r['fam'] != 'EC r148']), ('138 EC', [r for r in rows if r['fam'] == 'EC r148']), ('all 158', rows),
                      ('20 excl. genus 2', [r for r in rows if r['fam'] not in ('EC r148', 'GL4 genus 2')])):
    print('==', sel_name)
    base = None
    for name, m in models.items():
        s, n = chi2(m, sel)
        if 'no B' in name:
            s_law, n_law = chi2(lambda r: (-r['eps'] * r['c']) if r['A1'] else None, sel)
            print('   %-26s chi2 = %8.1f on %d   (law on the same rows %.1f)' % (name, s, n, s_law))
        else:
            print('   %-26s chi2 = %8.1f on %d' % (name, s, n))
    # power law c^gamma with free scale
    best = None
    for gam in np.linspace(-0.5, 2.5, 301):
        num = den = 0.0
        for r in sel:
            x = (-r['eps'] * r['c'] ** gam * r['W'])
            num += (r['G'] * np.conj(x)).real / r['var']; den += abs(x) ** 2 / r['var']
        s_ = num / den
        ch = sum(abs(r['G'] - s_ * (-r['eps'] * r['c'] ** gam * r['W'])) ** 2 / r['var'] for r in sel)
        if best is None or ch < best[0]: best = (ch, gam, s_)
        if abs(gam - 1) < 1e-9: ch1 = (ch, s_)
    # profile CI for gamma: chi2 within +1 of min (each complex residual counts both parts)
    cis = []
    for gam in np.linspace(-0.5, 2.5, 301):
        num = den = 0.0
        for r in sel:
            x = (-r['eps'] * r['c'] ** gam * r['W']); num += (r['G'] * np.conj(x)).real / r['var']; den += abs(x) ** 2 / r['var']
        s_ = num / den
        ch = sum(abs(r['G'] - s_ * (-r['eps'] * r['c'] ** gam * r['W'])) ** 2 / r['var'] for r in sel)
        if ch <= best[0] + 1: cis.append(gam)
    print('   power law A = -eps s c^g: best g = %.2f (chi2 %.1f, scale %.3f), 1-sigma range g in [%.2f, %.2f]; at g = 1: chi2 %.1f, scale %.3f'
          % (best[1], best[0], best[2], min(cis), max(cis), ch1[0], ch1[1]))
    # pooled slope with se
    num = sum((r['G'] * np.conj(-r['eps'] * r['c'] * r['W'])).real for r in sel)
    den = sum(abs(r['c'] * r['W']) ** 2 for r in sel)
    se = math.sqrt(sum(abs(r['c'] * r['W']) ** 2 * r['var'] for r in sel)) / den
    numc = sum((r['Gc'] * np.conj(-r['eps'] * r['c'] * r['W'])).real for r in sel)
    print('   pooled slope (|pred|^2 weights) %.4f +- %.4f ; with the prime-square term removed (level 1 only) %.4f' % (num / den, se, numc / den))
    # inverse-variance weighted slope
    num = sum((r['G'] * np.conj(-r['eps'] * r['c'] * r['W'])).real / r['var'] for r in sel)
    den = sum(abs(r['c'] * r['W']) ** 2 / r['var'] for r in sel)
    numc = sum((r['Gc'] * np.conj(-r['eps'] * r['c'] * r['W'])).real / r['var'] for r in sel)
    print('   inverse-variance slope %.4f +- %.4f ; prime-square removed %.4f' % (num / den, 1 / math.sqrt(den), numc / den))
# by family, inverse-variance, with P2 removal
print('== by family (inverse-variance slope, se; P2-corrected)')
for fam in sorted(set(r['fam'] for r in rows)):
    sel = [r for r in rows if r['fam'] == fam]
    num = sum((r['G'] * np.conj(-r['eps'] * r['c'] * r['W'])).real / r['var'] for r in sel)
    numc = sum((r['Gc'] * np.conj(-r['eps'] * r['c'] * r['W'])).real / r['var'] for r in sel)
    den = sum(abs(r['c'] * r['W']) ** 2 / r['var'] for r in sel)
    # also |pred|^2 weighted (proposer's pooling) with se
    num2 = sum((r['G'] * np.conj(-r['eps'] * r['c'] * r['W'])).real for r in sel); den2 = sum(abs(r['c'] * r['W']) ** 2 for r in sel)
    se2 = math.sqrt(sum(abs(r['c'] * r['W']) ** 2 * r['var'] for r in sel)) / den2
    print('   %-24s n=%3d  inv-var %.4f +- %.4f (P2-corr %.4f) | |pred|^2-weighted %.4f +- %.4f' % (fam, len(sel), num / den, 1 / math.sqrt(den), numc / den, num2 / den2, se2))
# deficit families jointly: GL2 level 1 + GL3 Sym2 + EC r148
for name, fams in (('level1+Sym2+EC', ('GL2 level 1', 'GL3 Sym2 E', 'EC r148')), ('all degree>=2 primitive', ('GL2 level 1', 'GL3 Sym2 E', 'EC r148', 'GL2 Artin S3', 'GL4 genus 2'))):
    sel = [r for r in rows if r['fam'] in fams]
    for corr in (False, True):
        num = sum(((r['Gc'] if corr else r['G']) * np.conj(-r['eps'] * r['c'] * r['W'])).real / r['var'] for r in sel)
        den = sum(abs(r['c'] * r['W']) ** 2 / r['var'] for r in sel)
        print('   joint %-26s P2corr=%s: %.4f +- %.4f  (%.1f se below 1)' % (name, corr, num / den, 1 / math.sqrt(den), (1 - num / den) * math.sqrt(den)))
