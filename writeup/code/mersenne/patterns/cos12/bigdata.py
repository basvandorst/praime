"""Numbers behind the big pictures of the cosine rows, period 12, and the two charts that go with them.
1. The pictures (width w, numbers 1..N): which columns hold primes, how the lane colours repeat, the twin tracks.
2. The lane race to M = 10^7: D(x) = (primes p <= x, 5 <= p, p = 5 or 7 mod 12) - (p = 1 or 11 mod 12).
3. The Mersenne exponents from 5 on, by lane.
Writes bigdata.txt (printed) and the SVG charts race.svg and rails.svg.
usage: bigdata.py"""
import math, numpy as np
from scipy.special import expi
halfli = lambda x: 0.5 * (expi(0.5 * np.log(x)) - expi(math.log(2)))
from math import gcd

M = 10_000_000
s = np.ones(M + 3, bool); s[:2] = False
for i in range(2, int((M + 2) ** .5) + 1):
    if s[i]: s[i * i::i] = False
pr = np.flatnonzero(s[:M + 1]); pr = pr[pr >= 5]
up = (pr % 12 == 1) | (pr % 12 == 11)
out = []
def say(*a):
    t = ' '.join(str(x) for x in a); print(t); out.append(t)

say('== 1. The pictures')
def lane(r):
    return 'up' if r % 12 in (1, 11) else ('down' if r % 12 in (5, 7) else '-')
for w, N in ((240, 100000), (277, 100000), (100, 40000)):
    lit = [c for c in range(w) if gcd(c + 1, w) == 1]
    say(f'width {w} (= {w % 12} mod 12), numbers 1..{N:,}: {len(lit)} of {w} columns can hold primes beyond the first row')
    if w % 12 == 0:
        cols = ''.join({'up': 'U', 'down': 'D'}[lane(c + 1)] for c in lit)
        per = next(k for k in range(1, len(lit) + 1) if len(lit) % k == 0 and all((lit[(i + k) % len(lit)] - lit[i]) % w == (lit[k % len(lit)] - lit[0]) % w and lane(lit[(i + k) % len(lit)] + 1) == lane(lit[i] + 1) for i in range(len(lit))))
        say(f'  every column keeps one lane; colours of the {len(lit)} columns left to right: {cols}')
        say(f'  the pattern repeats every {lit[per % len(lit)] - lit[0]} columns ({w // (lit[per % len(lit)] - lit[0])} times across)')
        tr = [a for a in range(1, w + 1) if a % 12 in (5, 11) and gcd(a, w) == 1 and gcd(a + 2, w) == 1]
        say(f'  twin tracks (column pairs a, a + 2 that can both hold primes, a = 5 or 11 mod 12): {len(tr)},'
            f' {sum(a % 12 == 11 for a in tr)} up and {sum(a % 12 == 5 for a in tr)} down; {sum(a + 2 > w for a in tr)} of them wraps to the next row')
        p = pr[pr <= N]; c = (p - 1) % w; cnt = np.bincount(c, minlength=w)[lit]
        say(f'  primes 5..{N:,} per column: min {cnt.min()}, max {cnt.max()}, mean {cnt.mean():.1f}; up {int(up[pr <= N].sum())}, down {int((~up)[pr <= N].sum())}')
        tw = p[s[p + 2]]; a = (tw - 1) % w + 1; on = np.isin(a, tr)
        say(f'  twin pairs (p, p + 2) with 5 <= p <= {N - 2:,}: {len(tw)}; on a twin track: {int(on.sum())}; the others: {[(int(q), int(q) + 2) for q in tw[~on]]} (first row, 5 divides {w})')
    else:
        k = w % 12
        say(f'  row r + 1 starts {k} places later in the wave than row r, so a lane class moves {k if k <= 6 else k - 12:+d} columns per row (mod 12)')
        if w % 12 == 4:
            lanes = {(c, r % 3): lane(w * r + c + 1) for c in lit for r in range(3)}
            say(f'  column 0 over three rows: {[lanes[(0, r)] for r in range(3)]}; every column cycles with period 3 rows (100 = 4 mod 12, 3 * 4 = 12)')

say('\n== 2. The lane race to 10^7')
D = np.cumsum(np.where(up, -1, 1))
say(f'primes 5..10^7: {len(pr):,}; up {int(up.sum()):,}, down {int((~up).sum()):,}; D(10^7) = {int(D[-1])}')
lead = D < 0
say(f'the up lane leads (D < 0) after {int(lead.sum()):,} of {len(pr):,} primes ({100 * lead.mean():.2f}%)')
st = np.flatnonzero(lead[1:] & ~lead[:-1]) + 1
if lead[0]: st = np.r_[0, st]
en = np.flatnonzero(~lead[1:] & lead[:-1]) + 1
if len(en) < len(st): en = np.r_[en, len(pr)]
runs = [(int(pr[a]), int(pr[b - 1]), b - a, int(D[a:b].min())) for a, b in zip(st, en)]
say(f'stretches where the up lane leads: {len(runs)}')
big = sorted(runs, key=lambda r: -r[2])
for r in sorted(big[:8]): say(f'  from {r[0]:,} to {r[1]:,}: {r[2]:,} primes, deepest D = {r[3]}')
say(f'first lead at p = {runs[0][0]:,}; last lead below 10^7 ends at p = {runs[-1][1]:,}')
for x in (10 ** 4, 10 ** 5, 5 * 10 ** 5, 10 ** 6, 10 ** 7):
    i = np.searchsorted(pr, x, 'right') - 1
    say(f'  D({x:,}) = {int(D[i])}; sqrt(x)/ln(x) = {math.sqrt(x) / math.log(x):.1f}; Li(sqrt x)/2 = {float(halfli(x)):.1f}')
say(f'largest D below 10^7: {int(D.max())} at p = {int(pr[D.argmax()]):,}')
g = np.exp(np.linspace(math.log(1e3), math.log(M), 200001)); gi = np.searchsorted(pr, g, 'right') - 1
ratio = D[gi] * np.log(g) / np.sqrt(g); r2 = D[gi] / halfli(g)
say(f'D(x) / (sqrt(x)/ln x), averaged over ln x from 10^3 to 10^7: {ratio.mean():.3f}; range {ratio.min():.2f} .. {ratio.max():.2f}')
say(f'D(x) / (Li(sqrt x)/2), averaged the same way: {r2.mean():.3f} (prediction: 1); range {r2.min():.2f} .. {r2.max():.2f}')
for a_, b_ in ((1e3, 1e5), (1e5, 1e7)):
    k = (g >= a_) & (g < b_); say(f'  {a_:.0e} .. {b_:.0e}: mean {ratio[k].mean():.3f} and {r2[k].mean():.3f}')

say('\n== 3. The exponents by lane')
EXP = [5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
eu = [p for p in EXP if p % 12 in (1, 11)]; ed = [p for p in EXP if p % 12 in (5, 7)]
say(f'{len(EXP)} exponents from 5 to 136,279,841: up {len(eu)}, down {len(ed)}')
say(f'  up: {", ".join(f"{p:,}" for p in eu)}')
w_ = lambda p: math.exp(0.5772156649) * math.log2((6 if p % 4 == 1 else 2) * p) / p
say(f'  to no. 51 (82,589,933): up {sum(p % 12 in (1, 11) for p in EXP[:-1])} of {len(EXP) - 1}')

open('bigdata.txt', 'w').write('\n'.join(out) + '\n')

# chart 1: the race, min/max envelope of D per bin of 10^7/1000
Wd, Ht, L, R, T, B = 1000, 320, 56, 16, 20, 40
nb = 1000; edges = np.searchsorted(pr, np.linspace(0, M, nb + 1)[1:], 'right')
lo = np.minimum.reduceat(D, np.r_[0, edges[:-1]]); hi = np.maximum.reduceat(D, np.r_[0, edges[:-1]])
ymax = int(math.ceil(max(D.max(), 1) / 50) * 50); ymin = int(math.floor(min(D.min(), 0) / 50) * 50)
X = lambda x: L + (Wd - L - R) * x / M
Y = lambda d: T + (Ht - T - B) * (ymax - d) / (ymax - ymin)
xs = [(i + 0.5) * M / nb for i in range(nb)]
b0 = Y(0)
up_area = f'{X(xs[0]):.1f},{b0:.1f} ' + ' '.join(f'{X(x):.1f},{Y(max(h, 0)):.1f}' for x, h in zip(xs, hi)) + f' {X(xs[-1]):.1f},{b0:.1f}'
dn_area = f'{X(xs[0]):.1f},{b0:.1f} ' + ' '.join(f'{X(x):.1f},{Y(min(l, 0)):.1f}' for x, l in zip(xs, lo)) + f' {X(xs[-1]):.1f},{b0:.1f}'
ref = ' '.join(f'{X(x):.1f},{Y(float(halfli(x))):.1f}' for x in np.linspace(2e4, M, 300))
svg = [f'<svg viewBox="0 0 {Wd} {Ht}" role="img" aria-label="Lane race to ten million: down-lane primes minus up-lane primes" style="width:100%;height:auto;display:block">']
for g in range(ymin, ymax + 1, 50 if ymax - ymin <= 400 else 100):
    svg.append(f'<line x1="{L}" x2="{Wd - R}" y1="{Y(g):.1f}" y2="{Y(g):.1f}" stroke="{"#5b6b82" if g == 0 else "#1f2e44"}" stroke-width="1"/><text x="{L - 6}" y="{Y(g) + 4:.1f}" fill="#9fb0c6" font-size="12" text-anchor="end">{g}</text>')
for k in range(0, 11, 2):
    svg.append(f'<text x="{X(k * 10 ** 6):.1f}" y="{Ht - B + 18}" fill="#9fb0c6" font-size="12" text-anchor="middle">{k} M</text>')
svg.append(f'<polygon points="{up_area}" fill="#f59e0b" fill-opacity="0.8"><title>down lane ahead</title></polygon>')
svg.append(f'<polygon points="{dn_area}" fill="#38bdf8"><title>up lane ahead</title></polygon>')
for r in runs:
    if r[2] >= 40:
        svg.append(f'<line x1="{X(r[0]):.1f}" x2="{X(r[0]):.1f}" y1="{b0:.1f}" y2="{Y(r[3]):.1f}" stroke="#38bdf8" stroke-width="3"><title>up lane ahead from {r[0]:,} to {r[1]:,} ({r[2]:,} primes, deepest {r[3]})</title></line>')
svg.append(f'<polyline points="{ref}" fill="none" stroke="#edf4ff" stroke-width="2" stroke-dasharray="6 5"><title>average lead predicted by Chebyshev’s bias: Li(√x)/2, about √x/ln x</title></polyline>')
svg.append(f'<text x="{X(8.2e6):.1f}" y="{Y(float(halfli(8.2e6))) - 12:.1f}" fill="#edf4ff" font-size="12" text-anchor="middle">predicted average: ½ Li(√x)</text>')
fr = runs[0]
svg.append(f'<text x="{X(fr[0]) + 6:.1f}" y="{b0 + 30:.1f}" fill="#d6e2f2" font-size="12">up lane ahead (blue): first at {fr[0]:,}, {len(runs)} stretches to 10 million</text>')
svg.append(f'<text x="{Wd - R}" y="{T - 6}" fill="#9fb0c6" font-size="12" text-anchor="end">down-lane primes minus up-lane primes</text></svg>')
open('race.svg', 'w').write('\n'.join(svg))

# chart 2: the exponents on the two rails, log scale
Wd, Ht, L, R = 1000, 170, 70, 20
X = lambda p: L + (Wd - L - R) * (math.log10(p) - 0.5) / (8.5 - 0.5)
yu, yd = 50, 120
svg = [f'<svg viewBox="0 0 {Wd} {Ht}" role="img" aria-label="Mersenne exponents from 5 to 136,279,841 on the up and down lanes" style="width:100%;height:auto;display:block">']
for yv, name, col in ((yu, 'up', '#38bdf8'), (yd, 'down', '#f59e0b')):
    svg.append(f'<line x1="{L}" x2="{Wd - R}" y1="{yv}" y2="{yv}" stroke="{col}" stroke-opacity="0.45" stroke-width="2"/><text x="{L - 10}" y="{yv + 4}" fill="#d6e2f2" font-size="13" text-anchor="end">{name}</text>')
for e in range(1, 9):
    svg.append(f'<line x1="{X(10 ** e):.1f}" x2="{X(10 ** e):.1f}" y1="{yu - 22}" y2="{yd + 14}" stroke="#1f2e44"/><text x="{X(10 ** e):.1f}" y="{Ht - 8}" fill="#9fb0c6" font-size="12" text-anchor="middle">10^{e}</text>'.replace('10^', '10<tspan dy="-5" font-size="9">').replace(f'{e}</text>', f'{e}</tspan></text>'))
for i, p in enumerate(EXP):
    yv = yu if p % 12 in (1, 11) else yd; jit = (-8, 0, 8)[i % 3]
    svg.append(f'<circle cx="{X(p):.1f}" cy="{yv + jit}" r="5" fill="#f472b6" stroke="#0b1220" stroke-width="2"><title>{p:,} = {p % 12} mod 12 ({"up" if yv == yu else "down"})</title></circle>')
svg.append(f'<text x="{Wd - R}" y="16" fill="#9fb0c6" font-size="12" text-anchor="end">up {len(eu)}, down {len(ed)} (about half each expected)</text></svg>')
open('rails.svg', 'w').write('\n'.join(svg))
