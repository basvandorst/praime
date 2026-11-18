"""Second round of section 20: reads section20.md as written by section20.py and adds the sign lemma, the folding lemma,
the new ring runs, the folded wheels, the three-cell wheel, the Mersenne statement and the table of proved and fast wheels.
Every number comes from tournament.json or from the output file of the script named next to it."""
import json, re, os
T = json.load(open('tournament.json')); f = lambda n: format(n, ',')
G = 'writeup/code/mersenne/strategies/ten/gears/tour/'
L = lambda n: '[%s](%s%s)' % (n, G, n)
s = open('section20.md', encoding='utf-8').read()
def rep(old, new):
    global s
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new)
def before(anchor, text):
    global s
    assert s.count(anchor) == 1, anchor[:60]
    s = s.replace(anchor, text + anchor)

AS = open('agrawal_sign.txt').read()
a5 = re.search(r'r = 5: home with x - 1 on r cells: (\d+) \(of these not a base-2 pseudoprime: (\d+)\)', AS).groups()
assert 'lemma holds on every ring tried: True' in AS and 'disagreements with plain Python: 0' in AS
FU = open('fold_units.txt').read()
nfold = re.search(r'folds \(r, e\): (\d+)', FU).group(1)
silent = eval(re.search(r'folds that say nothing \(u = 1\): (\[.*\])', FU).group(1))
assert 'exceptions: 0' in FU
fib = re.search(r'liars of the 5-cell ring below 1e10: (\d+) \| of these also liars of the Fibonacci wheel \(n divides F\(n -\+ 1\)\): (\d+)', FU).groups()
assert fib[0] == fib[1]
FC = open('fold_check.txt').read()
lost = re.search(r'times folds: (\d+) \| liars that a fold would have thrown away: (\d+)', FC).groups(); assert lost[1] == '0'
cmpf = re.search(r'numbers up to 2\^64: (\d+) \| different: (\d+)', FC).groups(); assert cmpf[1] == '0'
pipe = re.search(r'same liars as the plain full ring on (\d+) of (\d+) rings; (\d+) liars', FC).groups(); assert pipe[0] == pipe[1]
R3 = open('ring3_check.txt').read()
r3a = re.search(r'\(numbers that fit\): (\d+) \| rings where the fast program gives a different list: (\d+)', R3).groups(); assert r3a[1] == '0'
r3b = re.search(r'numbers up to 2\^64\): (\d+) \| different: (\d+)', R3).groups(); assert r3b[1] == '0'

rep("that there is never a mistake is open: Agrawal's conjecture for the ring, an unsettled Frobenius test for the two cells)",
    "that there is never a mistake is open: for the ring it is Agrawal's conjecture with x + 1 in place of x − 1, the same statement when the number of cells is even; for the two cells it is an unsettled Frobenius test)")

lem = f'''3. The sign. Agrawal's congruence is (x − 1)^n ≡ x^n − 1 (mod n, x^r − 1); Pascal's wheel is (x + 1)^n ≡ x^n + 1. For odd n they are the same statement on a ring with an even number of cells: replace x by −x. On an odd number of cells r they differ, and the link is: n is home with x − 1 on 2r cells exactly when it is home with x − 1 on r cells and with x + 1 on r cells. Proof: x^(2r) − 1 = (x^r − 1)(x^r + 1), and the two factors differ by 2, a unit mod n; so the ring of 2r cells is the product of the rings mod x^r − 1 and mod x^r + 1, and x → −x carries the second onto the first and x − 1 onto −(x + 1). Checked on every odd number below 10^6 that is not prime, for r = 3, 5, 7, 9, 11 ({L('agrawal_sign.py')} → {L('agrawal_sign.txt')}). Lemma 1 has no counterpart for x − 1 on an odd ring, where the sum of the cells is 0 after the first turn: of the {a5[0]} numbers below 10^6 that are home with x − 1 on 5 cells, {a5[1]} are not base-2 pseudoprimes. So the list reaches Pascal's form on every ring, and Agrawal's form on even rings only.
4. Folding by symmetry. Let r be prime and H a group of residues mod r. Renumbering the cells k → hk (x → x^h) carries the ring onto itself and a home position onto a home position, so a number that is home satisfies (x^h + 1)^n ≡ x^(hn) + 1 for every h. Multiply over h in H: u^n ≡ u(x^n) for u = ∏ (x^h + 1). This u does not change under H, so apart from cell 0 it has one value on each of the e = (r − 1)/|H| classes of H, and so has each of its powers: the test runs on a wheel of e cells (the sums of x^k over a class, taken mod 1 + x + … + x^(r−1)). A mistake on r cells is a mistake on every folded wheel. Exact facts for every prime r below 100 and every prime power e dividing (r − 1)/2, {nfold} folds ({L('fold_units.py')} → {L('fold_units.txt')}): mod 1 + x + … + x^(r−1) the product of all r − 1 factors x^h + 1 is 1, so u is a unit; and u itself is 1, so that the fold says nothing, when 2 lies in H (then ∏ (x^(2h) − 1) = ∏ (x^h − 1), and x^h + 1 is the quotient of the factors), and for these r in no other case. That happens for (r, e) = {", ".join("(%d, %d)" % x for x in silent)}.

'''
before('**The contest (computed).**', lem)

new = ''
LOW = {'25': 5, '27': 9, '32': 16}
if T['ring3']:
    parts = []
    for r in sorted(T['ring3'], key=int):
        c = T['ring3'][r]
        parts.append(('%s cells: %s tested' % (r, f(c['tested']))) if r not in LOW else ('%s cells, the candidates that fit %d cells and not %s (the others are covered by lemma 2): %s tested' % (r, LOW[r], r, f(c['tested']))))
    fs = sum(c['first_step'] for c in T['ring3'].values())
    new += f'''- More rings with a faster program ({L('ring3.c')}; sums of products are kept unreduced and reduced once per cell). Each candidate first runs on a ring of half the size: for an odd r the cells k and −k taken together, where ((x + 1)²/x)^n ≡ (x^n + 1)²/x^n must hold, and for 32 cells the ring of 16 cells with a sign (x^16 = −1). Only what passes goes to the full ring. {"; ".join(parts)}. {f(fs)} passed the first step and {sum(c['liars'] for c in T['ring3'].values())} came home. The program gives the known mistakes of numbers that fit ({f(int(r3a[0]))} on 20 rings below 10^10) and the same cells as plain Python on {r3b[0]} rings with numbers up to 2^64 ({L('ring3_check.py')} → {L('ring3_check.txt')}). The second program {L('ring2.c')} repeated these rings on a sample of 29,743 lines of the list: nothing came home, and for the prime rings it tested the same numbers ({L('ring2_sample.txt')}).
'''
if T['folded']:
    newr = [r for r in sorted(T['folded'], key=int) if r not in T['census'] and r not in T['ring3']]
    oldr = [r for r in sorted(T['folded'], key=int) if r in T['census'] or r in T['ring3']]
    def one(r):
        c = T['folded'][r]
        t = '%s cells: %s tested (%s' % (r, f(c['tested']), ', '.join('folded onto %s cells: %s' % (e, f(q['tested'])) for e, q in sorted(c['folds'].items(), key=lambda x: int(x[0]))))
        if c['nofold']: t += ', no fold applies and the half-size ring was used: %s' % f(c['nofold'])
        return t + ')'
    onepass = ''
    if os.path.exists('fold37_show_1.err') and open('fold37_show_1.err').read().strip().endswith('done'):
        sh = re.findall(r'fold (\d+) (\d+) passed by (\d+)', open('fold37_show_1.out').read())
        if len(sh) == 1 and T['fold_pass_total'] == 1: onepass = ' (%s, on the fold of %s cells onto %s, where it moves the cells with order 2; it is not home on the full ring; ' % (f(int(sh[0][2])), sh[0][0], sh[0][1]) + L('fold37_show_1.out') + ')'
    new += f'''- Folded wheels first ({L('fold.c')}). For a prime r each candidate was put on the smallest folded wheel of lemma 4 that says something and for which n is not in H; only what passes goes to the full ring (a class with no such fold goes to the half-size ring first). {"; ".join(one(r) for r in newr)}. Candidates that passed a folded wheel: {T['fold_pass_total']}{onepass}; came home: {sum(T['folded'][r]['liars'] for r in T['folded'])}.{(" The rings of " + ", ".join(oldr) + " cells were run again this way, with the same numbers of candidates as before and no pass: a second method for those rings.") if oldr else ""} Folding loses nothing: of the {f(int(lost[0]))} pairs (mistake below 10^10 of a number that fits, fold) every one passes its fold; the folded wheels agree with plain Python on the full ring in {f(int(cmpf[0]))} cases with numbers up to 2^64; and the whole chain, run on numbers that fit, returns exactly the {pipe[2]} known mistakes on {pipe[1]} rings ({L('fold_check.py')} → {L('fold_check.txt')}).
'''
nxt = [x for x in [5, 7, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53] if x not in T['rings_done']][0]
if T['ring3'] or T['folded']:
    new += f'''- Together: below 2^64 no number comes home on a ring of up to {T['all_r_upto']} cells that it does not fit (the first ring not run has {nxt} cells). By lemma 3 the same holds for Agrawal's own congruence on every even number of cells up to {T['even_upto']}.
'''
before('- Without lemma 1 and without the list:', new)

extra = ''
FA = None
for tag in ('1e8', '1e7'):
    if os.path.exists('folded_alone_%s.json' % tag) and os.path.exists('folded_alone_%s.txt' % tag):
        FA = json.load(open('folded_alone_%s.json' % tag)); fatag = tag; break
if FA:
    rows = "\n".join("| %d | %d | %s | %s | %s | %s |" % (x['r'], x['e'], " or ".join(map(str, x['orders_tested'])), f(x['tested']), f(x['passed']), ", ".join(map(str, x['first'])) if x['first'] else "none") for x in FA['folds'])
    tot = sum(x['passed'] for x in FA['folds']); o2 = sum(x['by_order'].get('2', 0) for x in FA['folds'])
    t3 = sum(x['tested'] for x in FA['folds'] if min(x['orders_tested']) >= 3); assert all(x['passed'] == 0 for x in FA['folds'] if min(x['orders_tested']) >= 3) and FA['passes_with_order_3_or_more'] == 0
    extra += f'''**Folded wheels alone (computed).** Without the doubling wheel and without the list, the folded wheels were put on every odd number below {fatag.replace('1e', '10^')} that is not prime and does not fit the ring; each number went to the fold chosen as above ({L('folded_alone.py')} → {L('folded_alone_%s.txt' % fatag)}).

| Ring | Cells after folding | Order of n on the fold | Tested | Came home | First mistakes |
|---|---|---|---|---|---|
{rows}

All {f(tot)} mistakes are numbers whose square lies in H ({f(o2)} of {f(tot)}): n moves the cells of the folded wheel round with order 2. No number that moves the cells with order 3 or more came home: nothing on the folds where every tested number has such an order ({f(t3)} tests on the folds onto {', '.join(sorted(set(str(x['e']) for x in FA['folds'] if min(x['orders_tested']) >= 3), key=int))} cells), and nothing of that kind on the others. This is the fit rule one level down: n² ≡ 1 (mod r) says that n acts on the whole ring with order at most 2. The first mistakes of the fold of 5 cells, 323 and 377, are the first mistakes of the Fibonacci wheel, for the reason given next.

'''
extra += f'''**Pascal's wheel on 5 cells contains the Fibonacci wheel (proved).** For r = 5 and H = {{1, 4}} the fold is u = (x + 1)(x⁴ + 1) = 2 + x + x⁴, and x + x⁴ is 2 cos 72° = φ − 1 for the golden ratio φ, so u = φ + 1 = φ². A number without the factor 5 that is home on 5 cells therefore satisfies φ^(2n) ≡ φ² for n ≡ ±1 and φ^(2n) ≡ φ^(−2) for n ≡ ±2 (mod 5), in the integers extended by φ, mod n. With √5 F(k) = φ^k − (−1/φ)^k and k = n ∓ 1, which is even, √5 F(k) φ^k = φ^(2k) − 1 ≡ 0; φ is a unit and 5 has no factor in common with n, so n divides F(n − 1), respectively F(n + 1). Every mistake of the ring of 5 cells is a mistake of the Fibonacci wheel, and by lemma 1 of the doubling wheel: all {fib[0]} below 10^10 are ({L('fold_units.txt')}). For n ≡ ±2 (mod 5) such a number is what the prize of 620 dollars asks for: a base-2 pseudoprime that divides F(n + 1). So on the numbers that do not fit it, the ring of 5 cells can only err if the prize problem has a counterexample, and the run of that ring above tested the same {f(T['pair']['tested'])} numbers as the run of the pair.

'''
if os.path.exists('three_check.json') and os.path.exists('three_check.txt'):
    TC = json.load(open('three_check.json')); TT = open('three_check.txt').read()
    pr = re.search(r'three primes just below 2\^64: (\d+) \| came home: (\d+)', TT).groups(); assert pr[0] == pr[1]
    cp = re.search(r'full ring: (\d+) numbers \| different verdicts: (\d+)', TT).groups(); assert cp[1] == '0'
    res = []
    al = [k for k in TC if k.startswith('alone_')]
    if al:
        k = sorted(al, key=lambda x: float(x.split('_')[1]))[-1]; c = TC[k]
        res.append('alone, on every odd number below %s that is not prime (%s numbers): %d came home; %s were turned away because a wheel size divides them, and %d, the cubes, have no wheel' % (k.split('_')[1].replace('1e', '10^'), f(c['tested']), c['home'], f(c['away']), c['nowheel']))
    if 'list_0' in TC and 'list_1' in TC:
        tot = {x: TC['list_0'][x] + TC['list_1'][x] for x in ('tested', 'home', 'away', 'nowheel')}; assert tot['tested'] == T['list_count']
        res.append('on all %s candidates below 2^64, that is together with the doubling wheel: %d came home (%s turned away, %d with no wheel); the largest wheel size used was %d' % (f(tot['tested']), tot['home'], f(tot['away']), tot['nowheel'], max(int(r) for p in ('list_0', 'list_1') for r in TC[p]['use'])))
    extra += f'''**A three-cell wheel (that every prime comes home is proved; the rest is computed).** Folding gives a machine that is not a ring: a cubic cousin of the two-cell wheel. Take a prime r ≡ 1 (mod 3) for which 2 is not a cube mod r ({", ".join(map(str, TC['sizes'][:6]))}, …) and let H be the cubes mod r. The folded wheel has three cells, one for each class of H. A turn multiplies by u, a fixed rule that replaces each cell by a sum of the three; for r = 7 it is multiplication by 2 + 2 cos(2π/7), a root of u³ − 5u² + 6u − 1. For n take the first such r that does not divide n and for which n is not a cube mod r; then n is home when after n turns the three cells show what they showed after the first turn, moved round by one or two places. So n acts with order 3: by the rule above this is a wheel that n does not fit. A number divisible by a wheel size tried on the way is turned away, and a cube has no wheel. Results ({L('three_check.py')} → {L('three_check.txt')}): {"; ".join(res)}. All {f(int(pr[0]))} primes tried came home, and the program agrees with plain Python on the full ring for {f(int(cp[0]))} numbers. In the terms of the literature this is a Frobenius test with a cyclic cubic polynomial.

'''
    TO = open('three_orders.txt').read(); TH = open('three_hunt.txt').read()
    o1 = re.search(r'primes below (\d+) that are not \+-1 mod 7: (\d+) \| ten smallest periods \(period, prime\): \[\((\d+), (\d+)\), \((\d+), (\d+)\), \((\d+), (\d+)\)', TO).groups(); assert 'primes whose period is below the prime itself: 0 ' in TO
    o2 = re.search(r'primes with a period below 20000: (\d+) \| products of two of them, or of three of the first 220, that meet the condition at every factor: (\d+)', TO).groups(); assert o2[1] == '0'
    sq = re.search(r'squares of the ([\d,]+) odd primes below 3,000,000 \(up to [\d,]+\): came home (\d+)', open('three_squares.txt').read()).groups()
    h1 = re.search(r'periods d up to (\d+) \(norms above 1e66 left out: d = \[\] \) \| primes not \+-1 mod 7 with such a period: (\d+)', TH).groups()
    extra += f'''An attempt to break it (computed; not a proof). A wheel with a unit is usually fooled by a product of primes on which the wheel has a short period: 323 = 17 · 19 fools the golden wheel that way. For a prime p that is not ±1 (mod 7) the three cells form one field mod p, and a number n with the factor p can only be home if the period d of u mod p divides n − p or n − p² (the first when n ≡ ±p, the second when n ≡ ±p² mod 7). For all {f(int(o1[1]))} such primes below 10^5 the period is larger than the prime (the smallest periods: {o1[2]} for p = {o1[3]}, {o1[4]} for p = {o1[5]}, {o1[6]} for p = {o1[7]}). Then no product pq of two of them, p < q, can be home: the condition at q asks that its period divides q − p, which is smaller than q. Among the {o2[0]} primes with a period below 20,000 no product of two or three meets the condition at every factor ({L('three_orders.py')} → {L('three_orders.txt')}). Primes with a period d must divide the norm of u^d − 1; PARI factors these norms for every d up to {h1[0]} and finds {h1[1]} such primes, 3 and 5 ({L('three_hunt.py')} → {L('three_hunt.txt')}). So the usual construction finds nothing here. A square p² is home mod p for every such prime, so it needs a separate look: the squares of the {sq[0]} odd primes below 3 · 10^6 were tried and {sq[1]} came home ({L('three_squares.txt')}). Prime factors that are ±1 (mod 7) were not treated this way; for them only the runs above speak.

'''
before('The two-cell wheel. Its step is the smallest odd prime c', extra)

i = s.index('The statement "home on r cells implies prime or n² ≡ 1 (mod r)" is Agrawal\'s conjecture.'); j = s.index('It is a computation, not a proof, and it rests on the completeness of the published list.', i) + len('It is a computation, not a proof, and it rests on the completeness of the published list.')
s = s[:i] + f'''Agrawal's conjecture says: if (x − 1)^n ≡ x^n − 1 (mod n, x^r − 1), then n is prime or n² ≡ 1 (mod r). Pascal's wheel is the same statement with x + 1, and by lemma 3 the two agree when the number of cells is even. The checks of Agrawal's form that we found reported: r < 100 with n < 10^10, r = 5 with n < 10^11, and a distributed search over 10^10 < n < 10^17. The runs above reach 2^64 ≈ 1.8 · 10^19 for Pascal's form on every ring up to {T['all_r_upto']} cells, and with it for Agrawal's form on every even number of cells up to {T['even_upto']}. For Agrawal's form on an odd number of cells they say nothing new, because a number that fools it need not be on the list. All of this is computation, not proof, and it rests on the completeness of the published list.''' + s[j:]

DM = open('duo_mersenne.txt').read()
m1 = re.search(r'odd prime exponents p below (\d+) : (\d+) \| steps used by the two-cell wheel: \[3\] \| two-cell wheel differs from Lucas-Lehmer: 0', DM).groups()
home = eval(re.search(r'exponents whose wheel came home: (\[.*\])', DM).group(1))
m2 = re.search(r"Pascal's wheel, the (\d+) exponents below (\d+) : differs from Lucas-Lehmer: 0", DM).groups()
bm = re.findall(r'p = (\d+) \((\d+) digits\): Lucas-Lehmer True ([\d.]+) \| two-cell wheel True ([\d.]+) \| Pascal on (\d+) cells True ([\d.]+) \| three-cell wheel from (\d+) cells True ([\d.]+) \| two cells with check turns True ([\d.]+)', DM)
b0 = [x for x in bm if x[0] == '4423'][0]; b1 = bm[-1]
PW = open('proved_wheels.txt').read()
w1 = re.search(r"Wilson's wheel, n = 2 \.\. (\d+) : wrong answers 0", PW).group(1)
w2 = re.search(r'the proved ring test, n = 2 \.\. (\d+) : wrong answers 0 \| primes that needed the rings: (\d+) \| largest machine: ring of (\d+) cells with (\d+) offsets for n = (\d+)', PW).groups()
w4 = re.search(r"Lucas's wheel, n = 2 \.\. (\d+) : wrong answers 0 \| not decided 0 \| largest step needed: (\d+) \(for n = (\d+)\)", PW).groups(); w5 = re.search(r'3 \* 2\^2208 \+ 1 \((\d+) digits\): proof found True with step (\d+)', PW).groups()
w6 = re.search(r'two cells with check turns, n = 2 \.\. (\d+) : wrong answers 0 \| not decided 0 \| largest start needed: (\d+) \(for n = (\d+)\)', PW).groups(); assert '2^521 - 1: proof found True with step 3 and start 1' in PW and '2^4423 - 1: proof found True with step 3 and start 1' in PW
CW = open('curve_wheel.txt').read(); cw = re.search(r'curve wheel, n = 2 \.\. (\d+) : wrong answers 0 \| not decided 0 \| primes proved (\d+)', CW).groups(); assert 'primes 53 (tick 241), Mersenne exponents 12 (tick 127), twin primes 33 (tick 241)' in CW
w3 = re.search(r'\(log2 n\)\^2 = ([\d,]+) cells and about sqrt\(r\) log2 n = ([\d,]+) offsets', PW).groups()
mers = f'''**Mersenne numbers (proved).** For n = 2^p − 1 with p an odd prime, the two-cell wheel is the Lucas–Lehmer test. Here n ≡ 1 (mod 3) and n ≡ 3 (mod 4), so (3/n) = −1 and the step is 3. Home means (1 + √3)^n ≡ 1 − √3, so (1 + √3)^(n+1) ≡ (1 − √3)(1 + √3) = −2. Now (1 + √3)² = 2(2 + √3) and n + 1 = 2^p, so (1 + √3)^(n+1) = 2^(2^(p−1)) (2 + √3)^(2^(p−1)), and 2^(2^(p−1)) ≡ 2 (mod n) because 2^p ≡ 1 (mod n) and 2^(p−1) ≡ 1 (mod p). So a wheel that is home gives (2 + √3)^(2^(p−1)) ≡ −1 (mod n). That is the Lucas–Lehmer condition, and it proves that n is prime: for a prime factor q of n the element 2 + √3 then has order exactly 2^p among fewer than q² elements, so q² > 2^p > n. On Mersenne numbers machine 11 cannot make a mistake. Check: for the {m1[1]} odd primes p below {f(int(m1[0]))} the wheel and the Lucas–Lehmer test agree, always with step 3; the wheel is home for p = {", ".join(map(str, home))} ({L('duo_mersenne.py')} → {L('duo_mersenne.txt')}). Pascal's wheel agrees with Lucas–Lehmer for the {m2[0]} exponents below {m2[1]}, but there is no such proof for it.

'''
before('**Cost.**', mers)
table = f'''**Proved or fast (computed; the two theorems are literature).** Timings are one run in Python with gmpy2 on 2^4423 − 1 (1,332 digits), from {L('duo_mersenne.txt')}; on 2^{b1[0]} − 1 ({f(int(b1[1]))} digits) Lucas–Lehmer, the two-cell wheel, the two cells with check turns, Pascal's wheel and the three-cell wheel took {b1[2]}, {b1[3]}, {b1[8]}, {b1[5]} and {b1[7]} seconds. The programs are plain and not tuned; only Lucas–Lehmer and Pascal's wheel use the cheap remainder of a Mersenne number.

| Wheel | Cells | For which n | Proved | Cost for 2^4423 − 1 |
|---|---|---|---|---|
| Wilson's wheel: one cell, multiplied by 1, 2, 3, … turn after turn; −1 after n − 1 turns | 1 | every n | yes | n − 1 turns, and none can be skipped |
| The proved ring test of Agrawal, Kayal and Saxena | more than {w3[0]} | every n | yes | about {w3[1]} offsets on that ring |
| Lucas's wheel: one cell, multiplied by a every turn; home after n − 1 turns and at no turn (n − 1)/q, q a prime factor of n − 1 | 1 | every n whose n − 1 can be factored | yes | cannot start: the prime factors of n − 1 = 2(2^4422 − 1) are needed |
| Two cells with check turns: the two-cell wheel from the start (a, 1), home after n turns and the round trip of n + 1 turns not finished at any turn (n + 1)/q | 2 | every n whose n + 1 can be factored, so every 2^p − 1 | yes | {b0[8]} s (step 3, start 1, one check turn) |
| Lucas–Lehmer | 1 | 2^p − 1 | yes | {b0[2]} s |
| Curve wheel (Goldwasser and Kilian): the points of a curve y² = x³ + ax + b mod n form a wheel whose number of positions m changes with a and b; one is picked whose m has a large prime factor q, and the wheel must be home after m turns and not after m/q | 1 point | every n | yes | machine 16 counts the positions one by one and stops at 3,000,000; programs that find m in a cleverer way (Atkin and Morain) reach far beyond: the largest proof we found reported has 109,297 digits |
| Two-cell wheel | 2 | every n (a square has no step) | only for 2^p − 1 | {b0[3]} s |
| Three-cell wheel | 3 | every n (a cube has no wheel) | no | {b0[7]} s (folded from {b0[6]} cells) |
| Pascal's wheel | 5 and more | every n | no | {b0[5]} s on {b0[4]} cells |

Wilson's wheel gave no wrong answer for n up to {f(int(w1))}, Lucas's wheel none for n up to {f(int(w4[0]))} (never undecided; largest step {w4[1]}, for n = {f(int(w4[2]))}) and it proves 3 · 2^2208 + 1 ({w5[0]} digits) with step {w5[1]}; the two cells with check turns gave none for n up to {f(int(w6[0]))} (never undecided; largest start {w6[1]}, for n = {f(int(w6[2]))}) and prove 2^521 − 1 and 2^4423 − 1 with step 3 and start 1; the curve wheel gave none for n up to {f(int(cw[0]))} and proved all {cw[1]} primes there ({L('curve_wheel.py')} → {L('curve_wheel.txt')}); and the proved ring test none for n up to {w2[0]} (its largest machine there: {w2[2]} cells and {w2[3]} offsets, for n = {w2[4]}) ({L('proved_wheels.py')} → {L('proved_wheels.txt')}). The wheel that is proved for every n with one cell is slow, the proved ring test is far too large, and the fast proofs by one wheel need something special: the form 2^p − 1, or the factors of n − 1 or n + 1. Elliptic curves remove that need for every n, at the price of a long search for a suitable curve and a chain of proofs. One small wheel with few turns, for every n and with a proof, is what is missing.

'''
before('**What this is and is not.**', table)
rep("the check of Pascal's wheel up to 2^64, and a small statement on prime powers", "the check of Pascal's wheel up to 2^64 on every ring up to %d cells, the folding of the ring onto small wheels with its link to the Fibonacci wheel and the prize problem, the remark that the two-cell wheel is the Lucas–Lehmer test on Mersenne numbers, and a small statement on prime powers" % T['all_r_upto'])
PC = open('page_check_10_15.txt').read(); assert len(re.findall(r'Lamps lit so far: 669\. Mistakes so far: 0\. script errors 0', PC)) == 7
R12 = open('rows_12_13.txt').read(); assert len(re.findall(r'primes 53 \(tick 241\), Mersenne exponents 12 \(tick 127\), twin primes 33 \(tick 241\); lamps of numbers that are not prime: 0', R12)) == 4
rep("The tape and the square jumps, which gave the same results as machine 1, were removed.", "Machines 12 to 16 are the three-cell wheel, Wilson's wheel, Lucas's wheel, the two cells with check turns and the curve wheel, built the same way; they give the same three results (" + L('rows_12_13.py') + " → " + L('rows_12_13.txt') + ", " + L('curve_wheel.txt') + "), and with 5,000 wheels the page's own code lights 669 lamps and makes no mistake on each of machines 10 to 16 (" + L('page_check_10_15.txt') + "). The test box of machine 14 proves 3 · 2^2208 + 1 and that of machine 15 proves 2^4423 − 1. The tape and the square jumps, which gave the same results as machine 1, were removed.")
open('section20.md', 'w', encoding='utf-8').write(s); print(len(s))
