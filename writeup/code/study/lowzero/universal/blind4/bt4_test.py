"""degree-4 follow-up: every registered reading with jackknife errors over L-functions and its decision.
usage: bt4_test.py         (predictions_blind4.json, zeros_blind4/)  -- refuses a prediction file whose sha256 != FROZEN
       bt4_test.py --dry   (dry/predictions_dry4.json, round-151 zero files) -- refuses unless sha256 == DRY
Excluded: keys without a COMPLETE zero file (listed). Writes results_blind4.json (results_dry4.json)."""
import sys, os, json, math, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from ulib import theta
FROZEN = 'e39d3eaaa4c54d438585c3c3a9c1e58581e3244ca5327c0be23853b6f4332467'
DRY = '1ddc51ee4754d7625c943151bb685a40ca0188874359df261fd3928978e79088'
R151Z = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'blind', 'zeros_blind')
dry = '--dry' in sys.argv
pf, zdir, want, outf = ('dry/predictions_dry4.json', R151Z, DRY, 'results_dry4.json') if dry else ('predictions_blind4.json', 'zeros_blind4', FROZEN, 'results_blind4.json')
raw = open(pf, 'rb').read(); h = hashlib.sha256(raw).hexdigest()
if h != want: sys.exit('REFUSED: sha256(%s) = %s, registered %s' % (pf, h, want))
P = json.loads(raw); print('prediction file %s sha256 %s (matches the registered hash)' % (pf, h))
LINES = {'2': 2.0, '1/2': 0.5, '3': 3.0, '1/3': 1 / 3}
per = {}; excluded = []
for key, e in P['L'].items():
    fn = os.path.join(zdir, key + '.npz')
    if not os.path.exists(fn): excluded.append([key, 'MISSING']); continue
    z = np.load(fn); info = json.loads(str(z['info']))
    if info.get('status') != 'COMPLETE' and not (dry and info.get('status') == 'REJECTED_ROUND151_RULE'):
        excluded.append([key, info.get('status')]); continue
    T2 = e['T2']; g = z['g']; g = g[(g > 0) & (g <= T2)]
    ph = np.exp(2j * theta(g, e['Q'], e['mu'])); w = np.sin(math.pi * g / T2) ** 2
    per[key] = dict(G=complex(ph.sum()), Gh=complex((w * ph).sum()), Gb=complex(ph[g > T2 / 4].sum()),
                    Gl={nm: complex((ph * np.exp(-1j * g * math.log(r))).sum()) for nm, r in LINES.items()}, n=len(g))
def xy(k, which):
    e = P['L'][k]; r = per[k]
    if which == 'r1': pairs = [(r['G'], complex(*e['pred']))]
    elif which == 'hann': pairs = [(r['Gh'], complex(*e['pred_hann']))]
    elif which == 'band': pairs = [(r['Gb'], complex(*e['pred_band']))]
    else: pairs = [(r['Gl'][nm], complex(*e['lines'][nm])) for nm in LINES]
    return sum((G * p.conjugate()).real for G, p in pairs), sum(abs(p) ** 2 for G, p in pairs)
def slope(keys, which='r1'):
    keys = [k for k in keys if k in per]
    if len(keys) < 3: return None
    X = np.array([xy(k, which) for k in keys]); sx, sy = X[:, 0].sum(), X[:, 1].sum()
    if sy == 0: return None
    b = sx / sy; loo = (sx - X[:, 0]) / (sy - X[:, 1]); n = len(keys)
    se = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    return dict(n=n, slope=b, se=se)
def dec(r): return None if r is None else ('PASS' if abs(r['slope'] - 1) <= 3 * r['se'] else 'FAIL')
OUT = dict(predictions=pf, sha256=h, excluded=excluded, readings={})
S = P['sets']
def show(name, r, extra=''):
    OUT['readings'][name] = dict(r, decision=dec(r)) if r else None
    print('%-34s %s %s' % (name, ('n=%3d  b = %.3f ± %.3f  %s' % (r['n'], r['slope'], r['se'], dec(r))) if r else 'no complete L', extra))
print()
show('R1 all degree 4', slope(S['all']))
rm, rp = slope(S['sym3_eps-1']), slope(S['sym3_eps+1'])
show('R2 Sym3 eps = -1', rm); show('R2 Sym3 eps = +1', rp)
if rm and rp:
    d = rp['slope'] - rm['slope']; sd = math.sqrt(rm['se'] ** 2 + rp['se'] ** 2)
    verdict = 'FLUCTUATION (law holds for eps = -1)' if (abs(rm['slope'] - 1) <= 3 * rm['se'] and rm['slope'] / rm['se'] >= 3) else \
              ('EFFECT' if (1 - rm['slope'] > 4 * rm['se'] and abs(rp['slope'] - 1) <= 3 * rp['se']) else 'UNDECIDED')
    OUT['readings']['R2 verdict'] = dict(difference=d, se=sd, verdict=verdict)
    print('   R2 difference b(+1) - b(-1) = %.3f ± %.3f ; watch-item verdict: %s' % (d, sd, verdict))
for s in ('S70', 'S50'):
    for nm in ('-1', '+1'):
        r = slope(S['sym3_%s_eps%s' % (s, nm)])
        OUT['readings']['R2s %s eps %s' % (s, nm)] = r
        print('   R2s %s eps = %s: %s' % (s, nm, ('n=%d b = %.3f ± %.3f' % (r['n'], r['slope'], r['se'])) if r else 'n < 3'))
show('R3 genus 2', slope(S['g2']))
for sg, nm in ((-1, '-1'), (1, '+1')):          # R3s: genus 2 by root number (reported, no decision; parallels the Sym^3 watch item)
    r = slope([k for k in S['g2'] if P['L'][k]['eps'][0] * sg > 0]); OUT['readings']['R3s genus 2 eps ' + nm] = r
    print('   R3s genus 2 eps = %s: %s' % (nm, ('n=%d b = %.3f ± %.3f' % (r['n'], r['slope'], r['se'])) if r else 'n < 3'))
show('R4 Rankin-Selberg', slope(S['rs']))
# R5 amplitude follows c (weighted LS, jackknife over L)
keys = [k for k in S['all'] if k in per]
def fit(ks):
    X = []; Y = []; W = []
    for k in ks:
        e = P['L'][k]; eps = complex(*e['eps']); unit = e['T2'] / (2 * math.pi)
        Y.append((per[k]['G'] * (-eps).conjugate()).real / abs(eps) / unit); X.append([1.0, e['c'] - 1.0]); W.append((unit / e['sd']['along']) ** 2)
    X = np.array(X); Y = np.array(Y); W = np.array(W)
    return np.linalg.solve(X.T @ (W[:, None] * X), X.T @ (W * Y))
if len(keys) >= 5:
    a, b = fit(keys); loo = np.array([fit([k for k in keys if k != j]) for j in keys]); n = len(keys)
    sa, sb = [math.sqrt((n - 1) / n * np.sum((loo[:, i] - loo[:, i].mean()) ** 2)) for i in (0, 1)]
    OUT['readings']['R5 amplitude'] = dict(alpha=a, se_alpha=sa, beta=b, se_beta=sb, decision='PASS' if abs(b - 1) <= 3 * sb else 'FAIL',
                                           H_c1='rejected' if b / sb >= 3 else 'NOT rejected')
    print('%-34s alpha = %.3f ± %.3f, beta = %.3f ± %.3f  %s ; H-c1 %s' % ('R5 amplitude follows c', a, sa, b, sb, OUT['readings']['R5 amplitude']['decision'], OUT['readings']['R5 amplitude']['H_c1']))
show('R6 lines', slope(S['all'], 'lines')); show('R7 Hann', slope(S['all'], 'hann'))
for nm, ks in (('all', S['all']), ('Sym3 eps = -1', S['sym3_eps-1']), ('Sym3 eps = +1', S['sym3_eps+1'])):
    r = slope(ks, 'band'); OUT['readings']['R8 band ' + nm] = r
    print('   R8 band (T2/4, T2] %-14s %s' % (nm, ('n=%d b = %.3f ± %.3f' % (r['n'], r['slope'], r['se'])) if r else 'n < 3'))
r1 = OUT['readings'].get('R1 all degree 4')
fal = bool(r1) and abs(r1['slope'] - 1) > max(4 * r1['se'], 0.10)
OUT['falsified'] = bool(fal)
print('\nexcluded:', excluded if excluded else 'none')
print('FALSIFIER (all new degree-4 slope outside 1 ± max(4 se_jk, 0.10)):', 'TRIGGERED' if fal else 'not triggered')
print('strict decisions failing:', [k for k, v in OUT['readings'].items() if isinstance(v, dict) and v.get('decision') == 'FAIL'] or 'none')
def cj(x):
    if isinstance(x, complex): return [x.real, x.imag]
    if isinstance(x, dict): return {k: cj(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [cj(v) for v in x]
    if isinstance(x, np.floating): return float(x)
    return x
OUT['per_L'] = cj({k: dict(G=v['G'], Gh=v['Gh'], Gb=v['Gb'], Gl=v['Gl'], n=v['n']) for k, v in per.items()})
json.dump(cj(OUT), open(outf, 'w'), indent=1)
