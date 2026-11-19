"""Writes section 21 of MERSENNE.md (the gaps between Mersenne exponents) to ~/scratch/prime3-lab/tour/section21.md.
Every number comes from mexp_gaps.json (mexp_gaps.py) and is checked by mexp_gaps_check.py."""
import json, os, re
J = json.load(open(os.path.expanduser('~/scratch/expgaps/mexp_gaps.json')))
C = open(os.path.expanduser('~/scratch/expgaps/mexp_gaps_check.txt')).read()
assert 'recomputed: 0 differ' in C
f = lambda n: format(n, ',')
G = 'writeup/code/mersenne/strategies/ten/expgaps/'
L = lambda n: '[%s](%s%s)' % (n, G, n)
s5, s127 = J['sizes_from_5'], J['sizes_from_127']
mods = J['mods']
W = J['wheels']
odd = [w for w in W['rows'] if w['p'] < 0.05]
gaps = J['gaps']
pairs = mods['pairs']
exact = {int(m): float(e) for m, e in re.findall(r'm = (\d+): \d+ gaps are multiples of m; model A \(simulated\) [\d.]+, exactly \d+/phi\(\d+\) = ([\d.]+)', C)}
rows = ''.join('| %d | %d | %.2f | %.1f | %.1f | %.2f |\n' % (r['m'], r['real'], exact[r['m']], r['A'], r['B'], r['pB']) for r in mods['rows'])
pc = lambda x: ('%.2f' % x).replace('-', '−')
t = f'''## 21. Gaps between Mersenne exponents (computed; the size law is the conjecture of Lenstra, Pomerance and Wagstaff; the remainder law is elementary)

The question: do the gaps between consecutive Mersenne exponents follow a rule, in their size or in their remainders, that random primes would not follow? The test is run on both pages: two layouts on [prime2.html](writeup/code/mersenne/prime2.html) and machine 24 on [prime3.html](writeup/code/mersenne/prime3.html).

**The data.** The 52 known exponents. GIMPS has checked every exponent below 83,200,483 twice (fact 6301 of [FACTS.md](FACTS.md)), so the first 51, up to 82,589,933, are certainly consecutive; for the gap from 82,589,933 to 136,279,841 that is not yet certain, and it is left out of every count. The 50 certain gaps: {', '.join(f(g) for g in gaps[:50])}. All but the first (2 → 3) are even.

**Sizes (computed against a conjecture).** The heuristic of Lenstra, Pomerance and Wagstaff says that log₂ of the exponents behaves like a random scatter (a Poisson process) with e^γ ≈ 1.781 points per doubling. Then log₂(p_(n+1)/p_n) is exponential with mean e^(−γ) = {s5['model']:.4f}. Tested with 20,000 runs of the model for each statistic ({L('mexp_gaps.py')} → {L('mexp_gaps.txt')}):

| Exponents | Gaps | Mean log₂ ratio | p | Kolmogorov–Smirnov distance | p | Correlation of one ratio with the next | p |
|---|---|---|---|---|---|---|---|
| from 5 on | {s5['n']} | {s5['mean']:.3f} | {s5['pMean']:.2f} | {s5['D']:.3f} | {s5['pD']:.2f} | {pc(s5['rho'])} | {s5['pRho']:.2f} |
| from 127 on | {s127['n']} | {s127['mean']:.3f} | {s127['pMean']:.2f} | {s127['D']:.3f} | {s127['pD']:.2f} | {pc(s127['rho'])} | {s127['pRho']:.2f} |

The largest ratio is 2^{s5['largest']:.2f} (127 → 521), the smallest 2^{s5['smallest']:.4f} (9,689 → 9,941). Counted per doubling from 2^2 to 2^26 the exponents are {', '.join(str(c) for c in J['doublings']['counts'])}: mean {J['doublings']['mean']:.2f}, variance {J['doublings']['var']:.2f}, so variance/mean {J['doublings']['var'] / J['doublings']['mean']:.2f}; a Poisson count is this even or more in {100 * J['doublings']['p_low']:.0f} % of runs. The sizes fit the model; a large jump is not followed by a small one more often than chance.

**Remainders (elementary, then computed).** m divides p_(n+1) − p_n exactly when the two exponents leave the same remainder mod m. For a random large prime the remainder mod m is spread evenly over the φ(m) remainders prime to m, so two independent ones agree with probability 1/φ(m). The values of m were fixed before the run: 3, 4, 5, 7, 8, 9, 11, 12, 13, 16, 24. Counted over the {pairs} gaps between the exponents from 61 on, against two models with {f(mods['sims'])} runs each. In model A every exponent is replaced by a random prime of about its size (window ×/÷ 1.5 below 10,000, ±5 % above). Model B uses only random primes for which 2^q − 1 has no factor 2kq + 1 with k ≤ {mods['K']}, as a Mersenne exponent has none:

| m | Gaps that are multiples of m | {pairs}/φ(m) | Model A | Model B | p (model B, two-sided) |
|---|---|---|---|---|---|
{rows}
{L('mexp_gaps_check.py')} → {L('mexp_gaps_check.txt')} recomputes the gaps, the sizes and every count with sympy and plain integers, and finds model A within 0.6 of {pairs}/φ(m) for each m (0 differ). The smallest p is {min(r['pB'] for r in mods['rows']):.2f}; with eleven values of m that is what chance gives. The largest deviation is for m = 4 and 5. For m = 4 the exponents are more often ≡ 1 (mod 4) than random primes ({J['class4']['real']} of {J['class4']['n']} from 61 on; model B gives {100 * J['class4']['modelB']:.1f} %). This tilt is known: a factor 2p + 1 of 2^p − 1 is possible only for p ≡ 3 (mod 4), and trial division up to k = {mods['K']} captures only part of it.

**The gear version (machine 24 on prime3.html).** A tick is one gap; wheel m is home when m divides the gap. Over the wheels 2 to 60 and the same {pairs} gaps, {len(odd)} wheels are off at the 5 % level in model B: {', '.join('%d (%d home, %.1f expected)' % (w['m'], w['real'], w['B']) for w in odd)}. They are two pairs, m and 2m, so really two findings: few gaps are multiples of 15 and many are multiples of 21. A sequence drawn from model B itself gives {W['null_mean']:.1f} such wheels on average, and {len(odd)} or more in {100 * W['p_global']:.0f} % of runs. That is not unusual, and the wheels 15 and 21 were not chosen in advance.

**What the gaps cannot do.** If the model holds and no exponent below 142,351,031 is missing (every exponent below it has been tested once, fact 6301), the next exponent has median {J['next']['median'] / 1e6:.0f} million and lies below {J['next']['q90'] / 1e6:.0f} million with 90 % chance (heuristic). The gaps give no sharper prediction, because they look like those of random primes chosen at the right rate.

**On the pages.**
- prime2.html, "Gaps between Mersenne exponents: one row per gap, on a doubling scale" (`expgaprows`). Each row runs from one exponent to the next, drawn with the primes in between, so the length of a row is log₂ of the ratio.
- prime2.html, "Gaps between Mersenne exponents: remainders mod P" (`expgapmod`). Only the exponents, one row each, in the column of their remainder mod P; with Links: In order, a vertical link is a gap that is a multiple of P. The note counts that over all 51 certain exponents for any P and gives the binomial p against 1/φ(P).
- Both layouts are in {LG('lx.js') if False else '[lx.js](writeup/code/mersenne/strategies/ten/lx.js)'} (added by {L('edit32.py')}).
- Screenshots: {L('shots/eg_rows_1e6.jpg')}, {L('shots/eg_mod4.jpg')}, {L('shots/eg_mod5.jpg')}, {L('shots/eg_mod12.jpg')}.
- prime3.html, machine 24: the status line counts the multiples of 3, 4, 5 and 7 among the gaps so far. {L('page_check_24.txt')} shows that at tick 50 it gives the counts of the table (21, 26, 6, 10).

What this is and is not: the size test repeats, with today's 51 exponents, a check that Wagstaff and others made with fewer. The remainder test is elementary. Neither finds anything beyond the model.

'''
open(os.path.expanduser('~/scratch/prime3-lab/tour/section21.md'), 'w', encoding='utf-8').write(t)
print(len(t))
