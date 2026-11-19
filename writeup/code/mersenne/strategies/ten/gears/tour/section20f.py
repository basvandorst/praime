"""Sixth round of section 20: adds the four machines that find primes before their tick (machines 20 to 23) to section20.md.
Every number comes from early_sum.json (early_summary.py) and inherit_sum.json (inherit_summary.py)."""
import json, math
E = json.load(open('early_sum.json'))
I = json.load(open('inherit_sum.json'))
f = lambda n: format(n, ',')
G = 'writeup/code/mersenne/strategies/ten/gears/tour/'
L = lambda n: '[%s](%s%s)' % (n, G, n)
LG = lambda n: '[%s](writeup/code/mersenne/strategies/ten/gears/%s)' % (n, n)
s = open('section20.md', encoding='utf-8').read()


def before(anchor, text):
    global s
    assert s.count(anchor) == 1, anchor[:60]
    s = s.replace(anchor, text + anchor)


def rep(old, new):
    global s
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new)


ee = E['early']
sf = ee['self']
dy = {t[0]: t[2] for t in E['dyadic']}
wk = E['work']
cnt_rows = ee['count']
ic = I['c']
icm = I['c_minus']
firsts = [x[2] for x in ic['cnt']]
law = I['law']
fit = I['fit']
py = I['python']
both = py['n - 1 and n + 1']
minus = py['n - 1 only']
PC = open('page_check_20_23.txt').read()
for need in ('machine 20 at tick 12: Lamps lit so far: 669. Mistakes so far: 0. script errors 0', 'machine 21 at tick 19: Lamps lit so far: 669. Mistakes so far: 0. script errors 0',
             'machine 22 at tick 5,000: Lamps lit so far: 669. Mistakes so far: 0. script errors 0', 'machine 23 at round 3: Lamps lit so far: 669. Mistakes so far: 0. script errors 0'):
    assert need in PC, need
assert firsts[1:] == [3, 19, 173, 21767, 1688173, 2223990533]


def sfr(T):
    return sf['rows'][T - 1]


tab = '| Machine, and what one tick is | Every prime below … is on after tick 3 | after tick 5 | after tick 10 | after tick 30 | All 669 primes below 5,000 on at tick |\n|---|---|---|---|---|---|\n'
tab += '| 10 to 19: wheel n makes one turn | 4 | 6 | 11 | 31 | 4,999 |\n'
tab += '| 21: one new prime gear, taken from the lamps | %s | %s | %s | %s | %d |\n' % (f(sfr(3)[3]), f(sfr(5)[3]), f(sfr(10)[3]), f(sfr(30)[3]), sf['all_5000'][0])
tab += '| 20: every wheel squares its cells | 16 | 64 | 2,048 | 2,147,483,648 | 12 |\n'
tab += '| 23: a round in which every dark wheel asks the lit ones | %s | %s | not known | not known | %d |\n' % (f(firsts[4]), f(firsts[6]), I['upto5000'][-1][0])
def nextprime(T):
    n = T + 1
    while any(n % d == 0 for d in range(2, int(n ** 0.5) + 1)):
        n += 1
    return n


cq = [nextprime(T) ** 2 for T in (3, 5, 10, 30)]
tab += '| 22: a gear changes counters; it counts and does not find | %s (counted) | %s (counted) | %s (counted) | %s (counted) | 4,999 for the lamps |\n' % tuple(f(x) for x in cq)
sq_tab = '| Tick T | 3 | 5 | 10 | 12 | 15 | 20 | 25 | 30 | 31 |\n|---|---|---|---|---|---|---|---|---|---|\n| Lamps (primes below 2^(T+1)) | ' + ' | '.join(f(dy[T]) for T in (3, 5, 10, 12, 15, 20, 25, 30, 31)) + ' |\n'
sf_tab = '| Tick | Gear put in | Every number decided below | New lamps | Lamps |\n|---|---|---|---|---|\n' + ''.join('| %d | %s | %s | %d | %s |\n' % (r[0], f(r[1]), f(r[3]), r[4], f(r[5])) for r in [sfr(1), sfr(2), sfr(3), sfr(5), sfr(10), sfr(19), sfr(30)] + [sf['tick_100']])
cn_tab = '| x | Counters | Prime gears | Moves | Primes up to x | Primes per move |\n|---|---|---|---|---|---|\n' + ''.join('| %s | %s | %s | %s | %s | %.2f |\n' % (f(x), f(n), f(g), f(m), f(p), p / m) for x, n, g, m, p in cnt_rows)
in_tab = '| Rounds | Primes up to %s | Smallest such prime |\n|---|---|---|\n' % f(ic['X']) + ''.join('| %d | %s | %s |\n' % (d, f(n), f(first)) for d, n, first in ic['cnt'][1:])
big_tab = '| Size | Primes | Mean rounds | Most | log₄(ln n) + %.2f |\n|---|---|---|---|---|\n' % law['both'][2]
cls = {a: (b, c) for a, b, c in ic['mean']}
for a in (10, 15, 20, 25, 30):
    big_tab += '| 2^%d to 2^%d (every prime) | %s | %.2f | – | %.2f |\n' % (a, a + 1, f(cls[a][0]), cls[a][1], math.log(math.log(1.5 * 2 ** a)) / math.log(4) + law['both'][2])
for b in I['big']:
    if b['digits'] >= 12:
        big_tab += '| %d digits (random) | %s | %.2f | %d | %.2f |\n' % (b['digits'], f(b['count']), b['mean'], b['worst'][0], math.log(math.log(5.5 * 10 ** (b['digits'] - 1))) / math.log(4) + law['both'][2])
m28 = dict((a, c) for a, b, c in ic['mean'])[28]
m28m = dict((a, c) for a, b, c in icm['mean'])[28]
h2 = lambda x: '%.2f' % (int(x * 100 + 0.5 + 1e-9) / 100)
box = E['box']
UNSEEN = ''
if I.get('unseen'):
    MORE = ''
    if I.get('more'):
        mm = I['more']
        assert max(mm, key=lambda m: m['mean'] - m['pred'])['digits'] == 13 and min(mm, key=lambda m: m['mean'] - m['pred'])['digits'] == 28
        assert min(((c - math.log(math.log(1.5 * 2 ** a_)) / math.log(4)), a_) for a_, b_, c in ic['mean'] if a_ >= 10)[1] == 17
        b40 = [b for b in I['big'] if b['digits'] == 40][0]
        assert b40['mean'] - math.log(math.log(5.5e39)) / math.log(4) > law['both'][2] + 0.02
        MORE = ('Because the large sizes allow only small samples, the rule was also tried on %d sizes between those of the fit, with %s random primes in all (' % (len(mm), f(sum(m['count'] for m in mm))) + L('inherit_big_more.out') + '):\n\n| Size | Primes | Mean rounds | Standard error | The rule says | Difference | In standard errors |\n|---|---|---|---|---|---|---|\n'
                + ''.join('| %d digits | %s | %.3f | %.3f | %.3f | %s | %s |\n' % (m['digits'], f(m['count']), m['mean'], m['se'], m['pred'], ('%+.3f' % (m['mean'] - m['pred'])).replace('-', '−'), ('%.1f' % m['z']).replace('-', '−')) for m in mm)
                + '\nHere the answer is sharp. The rule is never off by more than %.2f rounds, but it is not exact: %d of the %d differences are more than three standard errors. The mean moves in a slow wave around log₄(ln N) + %.2f, above it from 11 to 18 digits and below it near 28. So the honest form of the rule is: rounds = log₄(ln N) + %.2f, give or take %.2f. The crests of the wave (13 digits, and again 40 to 45 digits) and its troughs (the class 2^17, and again 28 digits) lie about one round apart on the scale of ln ln N, where a round is ln 4 = 1.39: that is what a count in whole rounds would do. With two crests and two troughs this is an impression, not a result.\n\n' % (max(abs(m['mean'] - m['pred']) for m in mm), sum(1 for m in mm if abs(m['z']) > 3), len(mm), law['both'][2], law['both'][2], max(abs(m['mean'] - m['pred']) for m in mm)))
    NEWC = ((I.get('second8') or {}).get('third') or {}).get('classes', [])
    DRIFT = ((I.get('second8') or {}).get('third') or {}).get('drift', [])
    zs = lambda z: ('%.1f' % z).replace('-', '−')
    rows = ''.join('| 2^%d to 2^%d, every prime | %s | %.3f | – | %.3f | – |\n' % (e['a'], e['a'] + 1, f(e['primes']), e['mean'], e['pred']) for e in NEWC)
    singles = []
    twice = None
    t75 = ([u for u in I['unseen'] if u['digits'] == 75 and len(u['parts']) > 1] or [None])[0]
    nall = all(len(u['parts']) > 1 for u in I['unseen'])
    for u in I['unseen']:
        if len(u['parts']) > 1:
            twice = u
            for k_, p_ in enumerate(u['parts']):
                rows += '| %d digits, random, %s sample | %d | %.3f | %.3f | %.3f | %s |\n' % (u['digits'], ('first', 'second', 'third')[k_], p_['count'], p_['mean'], p_['se'], u['pred'], zs(p_['z']))
                singles.append(abs(p_['mean'] - u['pred']))
            rows += '| %d digits, the samples together | %d | %.3f | %.3f | %.3f | %s |\n' % (u['digits'], u['count'], u['mean'], u['se'], u['pred'], zs(u['z']))
        else:
            rows += '| %d digits, random | %d | %.3f | %.3f | %.3f | %s |\n' % (u['digits'], u['count'], u['mean'], u['se'], u['pred'], zs(u['z']))
            singles.append(abs(u['mean'] - u['pred']))
    UNSEEN = ('The rule tried on sizes it was not made from (computed; ' + L('inherit_runs.sh') + ' → ' + ', '.join(L(x) for u in I['unseen'] for x in u['files']) + '). The constant %.2f comes from the classes up to 2^31 and the random sizes up to 60 digits. Afterwards:\n\n| Size | Primes | Mean rounds | Standard error | The rule says | Difference in standard errors |\n|---|---|---|---|---|---|\n' % law['both'][2]
              + rows
              + '\nNone of the random primes needed more than %d rounds.' % max(u['most'] for u in I['unseen'])
              + ((' The second samples need a word. The first sample of 75 digits lay %s standard errors from the rule, so a second one was drawn for that size with another seed; it lies %s standard errors on the other side. A second sample drawn because the first stood out is not a clean test%s' % (zs(t75['parts'][0]['z']), zs(t75['parts'][1]['z']), ('. Afterwards a second sample was drawn for every other size as well, so that no size is treated differently; all are shown.' if nall else ', so both are shown.'))) if t75 else '')
              + ' Taken one by one, the samples of random primes differ from the rule by at most %.2f rounds.' % max(singles)
              + ((' The two whole classes come from the later scan (' + L('inherit_third_17e9.out') + '); the lower one was partly, up to 4.2·10^9, in the first version of the fit. They lie %s rounds above the rule, and that is not noise: the difference mean − log₄(ln N) rises slowly through the classes 2^%d to 2^%d (%s), and is still rising at the top of the scanned range. Over all classes it stays between %.2f and %.2f; the rule is a rule of thumb with that much play, not a law with a fixed constant.' % (' and '.join('%.2f' % (e['mean'] - e['pred']) for e in NEWC), DRIFT[0][0], DRIFT[-1][0] + 1, ', '.join('%.3f' % d[1] for d in DRIFT), min(law['both'][0], min(d[1] for d in DRIFT)), max(law['both'][1], max(d[1] for d in DRIFT)))) if NEWC else '')
              + '\n\n' + MORE)
SRC = ''
if len(I.get('sources', [])) >= 3:
    s5 = ([d for d in I['sources'] if d['k'] == 5] or [None])[0]
    names = {1: 'N − 1', 2: 'N − 1 and N + 1 (machine 23)', 3: 'N − 1, N + 1 and N² + 1', 5: 'N − 1, N + 1, N² + 1, N² + N + 1 and N² − N + 1'}
    SRC = ('The reason tested where it had not been fitted (computed; ' + L('inherit_sources.py') + ' → ' + ', '.join(L('inherit_sources_%d.out' % d['k']) for d in I['sources']) + '). The same reasoning says more: with k sources of helpers each source has to supply 1/(2k) of the digits, so the rounds should grow like log_(2k)(ln N), with slope 1/ln(2k) against ln ln N. The lit part of N² + 1 is such a source, and for k = 5 also those of N² + N + 1 and N² − N + 1 (the proof steps for these numbers, wheels of four, three and six cells, are literature and were not built here; only the rounds were counted, with helpers below N). The slope 1/ln(2k) contains no fitted number. Every prime in six classes from 2^10 to 2^21 and random primes of 9 to 30 digits:\n\n| Sources of helpers | k | Slope measured | 1/ln(2k) |\n|---|---|---|---|\n'
           + ''.join('| %s | %d | %.2f | %.2f |\n' % (names[d['k']], d['k'], d['slope'], d['predicted_slope']) for d in I['sources']) + '\nFor k = 2 this fourth program gives the same class means as the C program.' + ((' For k = 5 the straight line does not fit: slope %.2f where the reasoning says %.2f. The curve there is a step: almost flat from 2^10 to 2^21 (%.2f to %.2f rounds) and rising after that. A second measure, chosen after seeing the step and therefore worth less: from the first class to the last the mean rises by %s rounds for k = %s, where the reasoning says %s. So the size of the effect is right within about 10 %% in all four cases, and the straight line is a fair picture only when one round is short compared with the range of sizes.' % (s5['slope'], s5['predicted_slope'], s5['flat'][0], s5['flat'][1], ', '.join('%.2f' % d['rise'] for d in I['sources']), ', '.join(str(d['k']) for d in I['sources']), ', '.join('%.2f' % d['pred_rise'] for d in I['sources']))) if s5 else '') + '\n\n')
dp = I['deep']
wd = dp.get('wide')
WIDE = (' A wider first stage (' + L('inherit_runs.sh') + ' → ' + L('inherit_deep_wide.out') + ', ' + L('inherit_deep_wide2.out') + '; every even s up to the smaller of r and 3.9·10^13 / r: %s candidates, %s primes of 6 rounds or more) gave no smaller one below %s.' % (f(wd['cands']), f(wd['six']), f(wd['bound']))) if wd else ''
s2 = dp['stage2']
r0 = s2['rows'][0]
DEEP = ('Primes of 7 rounds exist (computed; a search by construction, ' + L('inherit_deep.py') + ' → ' + L('inherit_deep_first.out') + ', ' + L(s2['file']) + '). Each neighbour of a prime of g + 1 rounds has a prime factor of g rounds or more (proved: otherwise that neighbour is lit completely by round g − 1, and one whole neighbour is larger than √N by itself). So candidates were built on the deep primes already known. From the %s primes r of 5 or 6 rounds below 4.2·10^9 (' % f(dp['deep5']) + L('inherit_deep_4e9.txt') + '), q = s · r ± 1 with even s up to %s gave %s primes of 6 rounds or more; from these, N = s · q ± 1 below %s gave %d primes of 7 rounds. The smallest is %s = %d × %s %s 1; the box of the page proves it in 7 rounds as well.' % (f(dp['s1']), f(dp['six']), f(s2['bound']), s2['found'], f(r0[0]), r0[2], f(r0[3]), '+' if r0[4] == '+' else '−') + WIDE + ' This search finds examples. It does not show that this is the smallest prime of 7 rounds' + ((': that one lies between %s (see the lemma used twice, above; %s if only the list of both programs is used) and %s' % (f(I['second8']['no7']), f(I['second8']['no7_both']), f(r0[0]))) if I.get('second8') else '') + '. Whether the number of rounds is unbounded is not proved.\n\n')
s2p = I.get('second')
SECOND = (' A second C program, written separately from the description of the rule (' + L('inherit_second.c') + ' → ' + L('inherit_second_4e9.out') + '; primes and factors from one sieve, no Miller–Rabin test), gives the same number of primes for every round up to %s and the same list of %s primes of 5 and 6 rounds.' % (f(s2p['X']), f(s2p['deep']))) if s2p else ''
s8 = I.get('second8')
if s8:
    c4, c8 = s8['chain4'], s8['chain8']
    SECOND += ' Run on to %s (' % f(s8['X']) + L('inherit_second_8e9.out') + '), it finds %s primes, the number of the counting recursion, and still none of 7 rounds; the primes of 6 rounds are %s. That list was also made without a scan (' % (f(s8['primes']), ', '.join(f(x) for x in s8['six'])) + L('inherit_lemma.py') + ' → ' + L('inherit_lemma.txt') + '): by the lemma under "Primes of 7 rounds" below, such a prime is s · q + 1 with q one of the primes of 5 or more rounds below 4.2·10^9, and trying all %s such numbers gives the same %d primes. The lemma reaches further when it is used on both neighbours (' % (f(s8['lemma_cands']), len(s8['six'])) + L('inherit_lemma_chain.py') + ' → ' + L('inherit_lemma_chain_4e9.txt') + ', ' + L('inherit_lemma_chain_8e9.txt') + '). A prime q up to 4B of g + 1 rounds has a prime factor r of g or more rounds in q − 1 = s · r and one in q + 1 = s′ · r′. If both were above B, then s = s′ = 2 and r′ = r + 1, two primes next to each other. So q = s · r ± 1 with r up to B. From the list up to 4.2·10^9, on which both programs agree: %d primes of 6 rounds up to %s (%s candidates), and from these no prime of 7 rounds up to %s. From the %s primes of 5 or more rounds up to 8.4·10^9 that the second program lists (' % (c4['six'], f(c4['top6']), f(c4['cands']), f(c4['top7']), f(c8['list'])) + L('inherit_deep_8e9.txt') + '): %d primes of 6 rounds up to %s, and no prime of 7 rounds up to %s.' % (c8['six'], f(c8['top6']), f(c8['top7']))
    th = s8.get('third')
    if th:
        SECOND += ' A scan confirms the count made without one: the second program with a smaller table (' + L('inherit_third.c') + ' → ' + L('inherit_third_17e9.out') + ') went up to %s, found %s primes (the number of the counting recursion), none of 7 rounds, and as primes of 6 rounds exactly the %d of the lemma. Its list of %s primes of 5 or more rounds (' % (f(th['X']), f(th['primes']), len(th['six']), f(th['deep'])) + L('inherit_deep_17e9.txt') + ', ' + L('inherit_lemma_chain_17e9.txt') + ') gives %d primes of 6 rounds up to %s and no prime of 7 rounds up to %s.' % (th['chain']['six'], f(th['chain']['top6']), f(th['chain']['top7']))
text = f'''**Primes before their tick: four machines that do not wait (squaring, the sieve, Legendre's count and the proof steps of Pocklington and Morrison are literature; the runs are computed; the rule log₄(ln n) for the rounds of machine 23 is a fit to the runs, and the reason given for it is a heuristic).** On machines 10 to 19 wheel n decides at tick n: 30 ticks give the 10 primes up to 30. Four machines of the page break that rule, each in another way. The table says below which number every prime is on after 3, 5, 10 and 30 ticks, with as many wheels as the machine can use ({L('early.py')} → {L('early.txt')}, {L('early_summary.py')} → {L('early_sum.txt')}; {L('inherit.py')} → {L('inherit.txt')}, {L('inherit_summary.py')} → {L('inherit_sum.txt')}).

{tab}
Machine 20, squaring wheels (the verdict is that of machine 19; the count of ticks is proved). The two cells of wheel n hold a + √c, the start of machine 19. Write n in binary as 1 b₁ b₂ … b_k. At tick i the wheel multiplies its two cells with themselves and, when b_i = 1, makes one more plain turn (it multiplies by a + √c).

1. After tick i the cells hold (a + √c)^m with m = (1 b₁ … b_i) in binary (proved: squaring doubles m and the plain turn adds b_i). So after k = ⌊log₂ n⌋ ticks they hold (a + √c)^n, and tick T decides every wheel below 2^(T+1).
2. No wheel that has to make n turns is earlier (proved). Let a tick be one multiplication of two values that the wheel already holds and at most one plain turn. Before the first tick the wheel holds 1 turn, fewer than 2. If the values held after t − 1 ticks stand for fewer than 2^t turns, the new value stands for at most 2(2^t − 1) + 1 < 2^(t+1) turns. So after t ticks fewer than 2^(t+1) turns are made, n turns need t ≥ ⌊log₂ n⌋ ticks, and the wheel is done exactly then. (A test that needs only the power (n ± 1)/2 could be one tick earlier.)
3. Checked ({L('early.py')}): for every n from 2 to {f(ee['ladder']['top'])} the verdict equals that of the stepping wheel ({L('norm2.py')}, which takes the power from the lowest binary digit up) and is right, and the number of ticks is the number of binary digits less one; differences {ee['ladder']['bad']} and {ee['ladder']['badt']}. In C ({L('norm2.c')}, {L('early_dyadic.sh')} → {L('early_dyadic.out')}) every odd number below 2^32 went on the wheel. The numbers that came home, added up:

{sq_tab}
   Every one of the 31 totals equals the count of machine 22's recursion for 2^(T+1) − 1, which shares nothing with the wheel.
4. Below 2^64 no number that is not prime comes home, given the list (see the two-cell wheel above). So at tick 63 the {f(E['pi_2_64'])} primes below 2^64 would all be on (that count is literature).

The price: every wheel works at every tick. With 5,000 wheels the odd wheels make {f(wk['ladder_ticks'])} ticks together, {wk['ladder_ticks'] / wk['primes']:.1f} for each prime found, and a tick is one squaring of two cells and at most one plain turn. Square-and-multiply is as old as the subject; the only thing added here is the clock.

Machine 21, gears from its own lamps (proved; the sieve of Eratosthenes counted in gears). Tick 1 puts in gear 2. At every tick after that the machine takes the smallest lamp that is not yet a gear.

1. Let the gears be the primes p₁, …, p_j and let q be the smallest number above p_j that no gear divides. Then q is prime: a prime factor of q below q is either at most p_j, and then it is a gear, or it lies between p_j and q, and then it is itself a number above p_j that no gear divides, smaller than q. A number m below q² that is not prime has a prime factor at most √m < q, and that factor is a gear. So every number below q² that no gear divides is prime (proved). Tick j lights the primes below p_(j+1)².
2. The run ({L('early.py')}; frontier up to {f(sf['top'])}, {sf['ticks']} ticks): at every tick each lamp is a prime and every prime below the frontier is lit, by comparison with a separate sieve.

{sf_tab}
3. An open question inside the machine. The new lamps of tick j ≥ 2 are the primes between p_j² and p_(j+1)²: 5, 6, 15, 9, 22, 11, 27, 47, 16, … (OEIS A050216). Brocard's conjecture (1904) says there are always at least four; it is open. In the run the fewest were {sf['min_new'][0]}, at tick {sf['min_new'][1]}.

The price: it looks at every number below the frontier. With 5,000 wheels that is {f(wk['strikes'])} strikes and {f(wk['looked'])} numbers looked at, {(wk['strikes'] + wk['looked']) / wk['primes']:.1f} moves for each prime.

Machine 22, counting wheels (proved; Legendre's rule of 1808 in the form known as Lucy's recursion, as in [p9_lucy.c](writeup/code/p9_lucy.c)). Let S(v, p) be the number of m with 2 ≤ m ≤ v that are prime or have no prime factor up to p. Then S(v, 1) = v − 1.

1. For a prime p with p² ≤ v the numbers that leave the count are p · m with p ≤ m ≤ v/p and m without a prime factor below p. There are S(⌊v/p⌋, p − 1) − π(p − 1) of them, and π(p − 1) = S(p − 1, p − 1). So S(v, p) = S(v, p − 1) − (S(⌊v/p⌋, p − 1) − S(p − 1, p − 1)) (proved). For p² > v nothing leaves, so S(v, p) = π(v) as soon as p ≥ √v: a counter is exact once the tick has passed the square root of its number.
2. p is prime exactly when S(p, p − 1) is one more than S(p − 1, p − 1) (proved: among the numbers up to p, only p itself can be without a prime factor below p, so S(p, p − 1) = π(p), and likewise S(p − 1, p − 1) = π(p − 1)). The machine needs no other test.
3. The counters needed for x are those of v = ⌊x/k⌋, because ⌊⌊x/k⌋/p⌋ = ⌊x/(kp)⌋: about 2√x of them. A move is one change of a counter.

{cn_tab}
   For each x the gears whose counter rose are exactly the primes up to √x, and the counters equal π ({L('early.py')}, against sympy). With the gears up to N and x = N² the moves stay below the primes counted for every N from {ee['cross']} on (tried up to 2,999). The box of the page counted the {f(box[1])} primes up to 10^10 in {f(box[2])} moves ({LG('early_node.js')}, {L('early_js_check.py')} → {L('early_js_check.txt')}).

The price: it knows how many primes there are, not which. A machine that finds its primes needs at least one move for each. The fastest counts of this kind need about x^(2/3) steps (Meissel, Lehmer, Lagarias, Miller and Odlyzko, Deléglise and Rivat; fact 2020 of [FACTS.md](FACTS.md)).

Machine 23, wheels that light each other (the proof steps are literature; the count of rounds is computed; its growth is a heuristic). Round 0: only the prime 2 is lit. In every round a dark number N looks at the lit primes that divide N − 1 and at those that divide N + 1. Let F₁ be the part of N − 1 made of lit primes (with their full powers), F₂ the same for N + 1, and L = lcm(F₁, F₂).

1. If L² > N, the lit primes are enough to decide N (proved, from steps 1 and 2 of the chain of wheels above; the combination is that of Brillhart, Lehmer and Selfridge, 1975). Run the one-cell test for the primes of F₁ and the two-cell test, with one step c for all of them, for the primes of F₂. If N fails one of them it is not prime. If it passes, every prime factor r of N is 1 mod F₁ and, with one sign for the whole of F₂, 1 or −1 mod F₂. So r is 1 or x mod L, where x is the number below L that is 1 mod F₁ and −1 mod F₂. A composite N has a prime factor r ≤ √N < L; r cannot be 1 mod L, so r = x. Hence N is not prime when 1 < x < N and x divides N, and prime otherwise: one division finishes the proof.
2. Every prime is lit in the end, and nothing else (proved). For a prime N the witnesses exist: a number that is not an f-th power mod N for the one cell, an element of order N + 1 among the elements of norm 1 for the two cells. By induction the primes below N are lit at some round, and when all prime factors of N − 1 are lit, L ≥ N − 1 > √N.
3. rounds(N) is 1 more than the smallest g such that the primes lit by round g are enough. Round 1 is known without a run (proved from the rule): the primes with (2^m)² > N for the power 2^m in N − 1 or N + 1, that is the Proth primes k · 2^m + 1 with k < 2^m and their mirror images k · 2^m − 1 with k < 2^m: 3, 5, 7, 11, 13, 17, 23, 31, 41, 47, 79, 97, …

The runs.

- Python with the real proofs ({L('inherit.py')} → {L('inherit.txt')}): every odd number up to {f(both['top'])} was put through the tests with witnesses. The numbers that got a proof are exactly the primes. As a test of step 1, {f(py['composites'][0])} odd composites below 3,000 were given every prime factor of n − 1 and n + 1 as helper: {py['composites'][1]} passed.
- C by the rule ({L('inherit.c')} → {L('inherit_c_4e9.out')}; primes by a Miller–Rabin test with the bases 2 to 13, which is exact far beyond this range) up to {f(ic['X'])}: {f(ic['primes'])} primes, the number that the counting recursion of machine 22 gives (the table leaves out the prime 2, which is lit from the start).{SECOND}

{in_tab}
- From the top down with PARI ({L('inherit_big.py')} → {L('inherit_big_check.txt')}): the same rule, worked out from N downwards, agrees on the {f(I['topdown'])} odd primes up to 200,000 and gives 6 rounds for 2,223,990,533. Its neighbours are 4 × 555,997,633 and 2 × 3 × 370,665,089, and both large primes need 5 rounds.

- A third count ({L('inherit_hist.py')} → {L('inherit_hist.txt')}): the rule from the top down for every prime up to {f(I['hist']['top'])}; the numbers of primes per round equal those of the C program ({f(I['hist']['primes'])} odd primes). It also confirms the {I['hist']['listed']} primes of 5 and 6 rounds that the C run printed; the three primes of 6 rounds up to 4.2·10^9 are {', '.join(f(x) for x in I['hist']['six'])}.
- The box of the page proves a typed number of up to 30 digits and lists its helpers round by round. For 2,223,990,533 it shows the helpers {', '.join(f(x) for x in I['boxcheck']['chain'])}; here 555,997,633 = 2 × 277,998,817 − 1 and 277,998,817 = 2 × 138,999,409 − 1. On {f(I['boxcheck']['primes'])} primes ({LG('inherit_node.js')}, {L('inherit_box_check.py')} → {L('inherit_box_check.txt')}) its rounds equal those of PARI, except for {I['boxcheck']['partial']} primes of {I['boxcheck']['digits'][1]} digits where the page cannot take a neighbour apart; there it says so and gives a longer proof.

With N − 1 alone ({L('inherit_c_minus_1e9.out')}, up to {f(icm['X'])}) the same count needs up to {icm['cnt'][-1][0]} rounds; the smallest primes of each round are {', '.join(f(x[2]) for x in icm['cnt'][1:])}. From 23 on, all 12 terms equal those of OEIS A056637 (the least prime of each class n− in the classification of Erdős and Selfridge, where every prime factor of p − 1 counts, down to 2 and 3). OEIS A082449 (the smallest prime for which taking the greatest prime factor of p − 1 again and again needs n steps) has one other term, 14,619,833 where these lists have 14,920,303. For the primes between 2^28 and 2^29 the mean is {m28m:.2f} rounds with N − 1 alone and {m28:.2f} with both neighbours.

How the rounds grow ({L('inherit_big.sh')} → {L('inherit_big.out')}: random primes, factored with PARI).

{big_tab}
The rule of thumb (heuristic). Suppose the primes below y are lit after g rounds. The part of a random number that is made of primes below y has, on average, about as many digits as y (the sum of ln p/(p − 1) over p ≤ y is ln y + O(1)). N is lit in the next round when F₁ and F₂ together have more than half the digits of N, and each neighbour supplies about ln y: 2 ln y > ½ ln N, that is N < y⁴. So one more round raises the size that can be reached to about its fourth power, and rounds(N) ≈ log₄(ln N) + constant. With N − 1 alone the same count gives N < y², hence log₂(ln N). The fit (computed): a straight line through {fit['points']} size classes (every class from 2^{I['fit_classes'][0]} to 2^{I['fit_classes'][1] + 1} and the random sizes from 12 to 60 digits; the table above shows some of them; mean rounds against ln ln N) has slope {fit['slope']:.2f}, where 1/ln 4 = 0.72; and with the slope fixed at 1/ln 4 the difference rounds − log₄(ln N) lies between {law['both'][0]:.2f} and {law['both'][1]:.2f} for every one of these classes. With N − 1 alone and the slope fixed at 1/ln 2 the difference lies between {str(round(law['minus'][0], 2)).replace('-', '−')} and {law['minus'][1]:.2f} for the classes 2^10 to 2^28.

{UNSEEN}{SRC}{DEEP}What is known, and what we did not find. The two proof steps are Pocklington's (1914) and Morrison's (1975). Counting levels of such proofs is known in two forms: the classes of Erdős and Selfridge (all prime factors of p + 1, or of p − 1, down to 2 and 3; OEIS A005113 and A056637) and the Pratt tree (all prime factors of p − 1), whose height Ford, Konyagin and Luca (2010) bound for almost all primes and for which their model suggests about e · ln ln p. The count used here differs: both neighbours, and only as many helpers as it takes to pass √N. We did not find it in a search, and the local copy of the OEIS (31 October 2026) contains neither 2, 3, 19, 173, 21767, 1688173, 2223990533 nor the list of round 1 as one sequence (its two halves are there: the Proth primes, A080076, and the primes k · 2^m − 1 with k < 2^m, A112715). It is a small step from known things and may well be known.

What machine 23 is not. It is not a faster way to find primes: a round is heavy, because every dark wheel tries in every round, and a machine that may turn as much as it likes within one tick can let each wheel decide alone, as machine 20 does. For a large N the hard part is taking N − 1 and N + 1 apart. What the rounds measure is how deep a proof by these two steps has to go.

On the page ({L('page_check_20_23.txt')}; the page's own code in headless Chrome with 5,000 wheels): machine 20 has 669 lamps at tick 12, machine 21 at tick 19, machine 23 after round 3, machine 22 at tick 5,000 with {f(wk['count_primes'])} primes counted; no mistake on any of them. The page code of the four machines was compared with Python ({LG('early_node.js')}, {L('early_js_check.py')} → {L('early_js_check.txt')}: {E['js']} comparisons, 0 differ; {LG('inherit_node.js')}, {L('inherit_js_check.py')} → {L('inherit_js_check.txt')}: rounds, helpers, witnesses and the division for {f(I['js'])} odd numbers, 0 differ).

'''
before('**What this is and is not.**', text)
rep("the chain of wheels and the test with a receipt are literature and their runs are computed)", "the chain of wheels and the test with a receipt are literature and their runs are computed; the four machines that find primes before their tick are literature with computed runs, and the rule for the rounds of machine 23 is a heuristic)")
open('section20.md', 'w', encoding='utf-8').write(s)
print(len(s))
