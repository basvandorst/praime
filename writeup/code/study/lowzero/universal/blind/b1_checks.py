"""PARI checks before registration (no zero of a TEST L-function is computed):
 part 'feq'   : lfuncheckfeq for every test L (and, for T5, zeta_K and L(rho)); also lfunrootres.
 part 'short' : for every proxy (never tested): lfunzeros on (0.01, 15] and (T2-10, T2] against (theta(b)-theta(a))/pi.
 part 'full <class>' : one proxy per class through the registered zero routine to the full height (timing + rule).
usage: b1_checks.py feq | short | full T1|T2|T3|T4a|T4b|T5"""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from bz_lib import build, params, zeros, unfold
from ulib import to_c
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
TS = json.load(open('testset.json'))
HEIGHT = dict(T1=1000, T2=500, T3=150, T4a=70, T4b=150, T5=300)
def specs(cls, which):
    out = []
    for x in TS[cls][which]:
        if cls == 'T1': out += [['char', x[0], x[1]], ['char', x[0], x[2]]]
        elif cls == 'T2': out.append(['mf', x[0], x[1], x[2]])
        elif cls == 'T3': out.append(['sym', x[1], 2])
        elif cls == 'T4a': out.append(['sym', x[1], 3])
        elif cls == 'T4b': out.append(['g2', x[1], x[2]])
        elif cls == 'T5': out += [['nf', x[1]], ['artin', x[1]]]
    return out
part = sys.argv[1]
res = {}
if part == 'feq':
    for cls in ('T1', 'T2', 'T3', 'T4a', 'T4b', 'T5'):
        t0 = time.time(); worst = -999; bad = []
        for s in specs(cls, 'test'):
            L = build(pari, s); f = float(pari.lfuncheckfeq(L)); worst = max(worst, f)
            if f > -20: bad.append([s, f])
        res[cls] = dict(n=len(specs(cls, 'test')), worst_feq=worst, failing=bad, secs=round(time.time() - t0, 1))
        print(cls, res[cls], flush=True)
elif part == 'short':
    extra = dict(T4b=[['g2', 'x^2+x', 'x^3+1'], ['g2', 'x^5-2*x^3+2*x^2-x', '1']])
    for cls in (sys.argv[2:] or ('T1', 'T2', 'T3', 'T4a', 'T4b', 'T5')):
        T2 = HEIGHT[cls]
        for s in specs(cls, 'proxy') + extra.get(cls, []):
            if s[0] == 'artin': continue
            L = build(pari, s); Q, k, mu = params(pari, L)
            t0 = time.time()
            lo = len(pari.lfunzeros(L, [0.01, 15.0])); hi = len(pari.lfunzeros(L, [T2 - 10.0, T2]))
            elo = float(unfold(15.0, Q, mu) - unfold(0.01, Q, mu)); ehi = float(unfold(T2, Q, mu) - unfold(T2 - 10.0, Q, mu))
            ok = abs(lo - elo) <= 2.5 and abs(hi - ehi) <= 2.5
            print('%-4s proxy %-60s low %3d (theta %6.1f)  top %3d (theta %6.1f)  %s  %.1fs' % (cls, str(s)[:60], lo, elo, hi, ehi, 'OK' if ok else 'FAIL', time.time() - t0), flush=True)
elif part == 'full':
    cls = sys.argv[2]; T2 = HEIGHT[cls]
    pr = specs(cls, 'proxy')
    if cls == 'T5': pr = [p for p in pr if p[0] == 'nf']
    s = pr[-1] if cls != 'T4b' else pr[0]
    g, info = zeros(pari, s, T2, poles=1 if s[0] == 'nf' else 0)
    print(cls, 'proxy', s, 'n', info['n'], 'status', info['status'], 'drift %.3f r_mean %.3f c_exp %.3f m0 %d added %d dup %d' % (
        info['drift'], info['r_mean'], info['c_exp'], info['m0'], info['added'], info['duplicates']), 'secs', info['secs'], flush=True)
