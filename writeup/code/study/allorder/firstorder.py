"""First-order (pair-correlation) slope without truncation over the range [H1, H2] = [1e4, 1e6]:
   s_inf = LR(y=17) - (kappa/lam) * Sigma,   Sigma = sum_{p > 17} f_p / p,
   f_p = clip(ln(H2/(c p)) / ln(H2/H1), 0, 1)   (turn-on of the power of p; c fitted to tailx1.log fractions),
   kappa/lam per direction = p * (step of the comparator from the previous y to y = p) / f_p  (round 139 report, lam 24.16)
   or = (p * inc(p))/lam from the cross-stretch increments (tailx: 0.34 at lam 22.35, scaled by (22.35/lam)^1.9)."""
import json, math, numpy as np
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)
PR = primes(10 ** 6)
T = json.load(open('tailx_19-47-97-199-499.json'))['turnon']
H1, H2 = 1e4, 1e6
def frac(p, c): return min(1.0, max(0.0, math.log(H2 / (c * p)) / math.log(H2 / H1)))
cs = np.linspace(50, 1500, 2901)
err = [sum((frac(int(p), c) - v['frac']) ** 2 for p, v in T.items()) for c in cs]
c = float(cs[int(np.argmin(err))])
print('turn-on fractions measured:', {p: round(v['frac'], 2) for p, v in T.items()}, ' fitted c = %.0f' % c, ' model:', {p: round(frac(int(p), c), 2) for p in T})
Sig = {cc: sum(frac(int(p), cc) / p for p in PR if p > 17) for cc in (c / 2, c, 2 * c)}
print('Sigma = sum_{p>17} f_p/p: c/2 %.3f  c %.3f  2c %.3f' % (Sig[c / 2], Sig[c], Sig[2 * c]))
rows = [l.split() for l in open('fresh139_cells.txt') if l[0].isdigit()]
print('cell   measured slope      LR17      kappa/lam(steps)  first-order s_inf   ratio meas/LR  ratio s_inf/LR')
out = {}
for r_ in rows:
    cell = r_[0]; sl = float(r_[8]); e = float(r_[10]); lr = float(r_[12]); lr13 = float(r_[14]); lr11 = float(r_[16])
    if cell.endswith('_50'): continue
    k13 = abs(lr - lr13) * 17 / frac(17, c)
    kap = k13 if math.isnan(lr11) else 0.5 * (k13 + abs(lr13 - lr11) * 13 / frac(13, c))
    s_inf = lr - kap * Sig[c]
    out[cell] = dict(meas=[sl, e], lr=lr, kappa_over_lam=kap, s_inf=s_inf, s_inf_range=[lr - kap * Sig[2 * c], lr - kap * Sig[c / 2]])
    print('%-5s  %+.5f +- %.5f  %+.5f   %.4f            %+.5f [%+.5f, %+.5f]   %.2f          %.2f' % (cell, sl, e, lr, kap, s_inf, lr - kap * Sig[2 * c], lr - kap * Sig[c / 2], sl / lr, s_inf / lr))
lamB = 22.35; kB = 0.34 / lamB
print('cross-stretch kappa/lam at lam %.2f: %.4f per direction; at 24.16 (x (22.35/24.16)^1.9): %.4f' % (lamB, kB, 0.34 * (22.35 / 24.16) ** 1.9 / 24.16))
json.dump(dict(c=c, Sigma={str(k): v for k, v in Sig.items()}, cells=out), open('firstorder.json', 'w'), indent=1)
