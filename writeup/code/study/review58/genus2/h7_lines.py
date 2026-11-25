"""Does the genus-2 deficit sit in the common constant c or only at r = 1?  Own line predictions for r = 2, 1/2, 3, 1/3:
G(r) = sum e^{2i theta(g)} r^{-ig} = -eps c (S_p(v)/F_p) T/2pi, S_p(v) = sum_{m - h = v} conj lambda(p^m) mu(p^h) p^{-(m+h)/2},
F_p = F_p(1/p) (fact 2055/2057 lines, as in the registered predictions), local factors at p = 2, 3 from lfuneuler (exact).
Checked against the registered line predictions of round 152 where they exist. Slopes pooled over the four lines (registered
estimator) for genus 2 (100), and the rotation Im-slope at r = 1."""
import sys, os, math, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from g2data import load, slope, B4
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
P4 = json.load(open(B4 + 'predictions_blind4.json'))['L']
rows = load()
def local(L, p):
    Pp = 1 / pari.lfuneuler(L, p); co = [float(c) for c in pari.Vecrev(Pp)] if str(pari.type(Pp)) == 't_POL' else [float(Pp)]
    co = [c / co[0] * p ** (-j / 2) for j, c in enumerate(co)]
    J = 60; lam = [1.0] + [0.0] * J
    for j in range(1, J + 1): lam[j] = -sum(co[i] * lam[j - i] for i in range(1, min(j, len(co) - 1) + 1))
    F = sum(lam[j] * co[j] * p ** (-j) for j in range(len(co)))
    S = {v: sum(lam[h + v] * co[h] * p ** (-(2 * h + v) / 2) for h in range(len(co)) if 0 <= h + v <= J) for v in (1, -1)}
    return F, S
LINES = {'2': (2, 1), '1/2': (2, -1), '3': (3, 1), '1/3': (3, -1)}
maxdiff = 0
for r in rows:
    L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2])); loc = {p: local(L, p) for p in (2, 3)}
    r['lp'] = {}; r['lg'] = {}
    for nm, (p, v) in LINES.items():
        F, S = loc[p]; r['lp'][nm] = -r['eps'] * r['c'] * S[v] / F * 150 / (2 * math.pi)
        rr = p ** v; r['lg'][nm] = complex((r['ph'] * np.exp(-1j * r['gg'] * math.log(rr))).sum())
        if r['key'] in P4: maxdiff = max(maxdiff, abs(r['lp'][nm] - P4[r['key']]['lines'][nm][0]))
print('own line predictions vs registered (round 152, 40 curves): max |diff| %.2e' % maxdiff)
for nm, S in (('all 100', rows), ('eps+1', [r for r in rows if r['eps'] > 0]), ('r151+152', [r for r in rows if r['rnd'] != 'rep']), ('rep', [r for r in rows if r['rnd'] == 'rep'])):
    Gs = [r['G'] for r in S]; Ps = [-r['eps'] * r['c'] * 150 / (2 * math.pi) for r in S]
    b1, s1 = slope(Gs, Ps)
    # pooled lines: per curve sum over the four lines (registered estimator pools pairs per L)
    x = np.array([sum((r['lg'][k] * np.conj(r['lp'][k])).real for k in LINES) for r in S]); y = np.array([sum(abs(r['lp'][k]) ** 2 for k in LINES) for r in S])
    n = len(S); bl = x.sum() / y.sum(); loo = (x.sum() - x) / (y.sum() - y); sl = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    per = []
    for k in LINES:
        bk, sk = slope([r['lg'][k] for r in S], [r['lp'][k] for r in S]); per.append('%s %.3f+-%.3f' % (k, bk, sk))
    rot = sum((g * p).imag for g, p in zip(Gs, Ps)) / sum(p * p for p in Ps)
    print('%-9s r=1 %.3f+-%.3f | lines pooled %.3f+-%.3f | %s | Im-slope at r = 1 %+.3f' % (nm, b1, s1, bl, sl, '  '.join(per), rot))
