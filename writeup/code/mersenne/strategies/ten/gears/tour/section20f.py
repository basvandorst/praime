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
    UNSEEN = 'The rule tried on sizes it was not made from (computed; ' + ', '.join(L('inherit_big_%d.out' % u['digits']) for u in I['unseen']) + '). The constant %.2f comes from the classes up to 60 digits. Afterwards: ' % law['both'][2] + '; '.join('%d random primes of %d digits needed %.3f rounds on average (standard error %.2f), and the rule says %.3f' % (u['count'], u['digits'], u['mean'], u['se'], u['pred']) for u in I['unseen']) + '. None needed more than %d rounds.\n\n' % max(u['most'] for u in I['unseen'])
SRC = ''
if len(I.get('sources', [])) >= 3:
    s5 = ([d for d in I['sources'] if d['k'] == 5] or [None])[0]
    names = {1: 'N − 1', 2: 'N − 1 and N + 1 (machine 23)', 3: 'N − 1, N + 1 and N² + 1', 5: 'N − 1, N + 1, N² + 1, N² + N + 1 and N² − N + 1'}
    SRC = ('The reason tested where it had not been fitted (computed; ' + L('inherit_sources.py') + ' → ' + ', '.join(L('inherit_sources_%d.out' % d['k']) for d in I['sources']) + '). The same reasoning says more: with k sources of helpers each source has to supply 1/(2k) of the digits, so the rounds should grow like log_(2k)(ln N), with slope 1/ln(2k) against ln ln N. The lit part of N² + 1 is such a source (its proof step, a wheel of four cells, is literature and was not built here; only the rounds were counted, with helpers below N). The slope 1/ln(2k) contains no fitted number. Every prime in six classes from 2^10 to 2^21 and random primes of 9 to 30 digits:\n\n| Sources of helpers | k | Slope measured | 1/ln(2k) |\n|---|---|---|---|\n'
           + ''.join('| %s | %d | %.2f | %.2f |\n' % (names[d['k']], d['k'], d['slope'], d['predicted_slope']) for d in I['sources']) + '\nFor k = 2 this fourth program gives the same class means as the C program.' + ((' For k = 5 the straight line does not fit: slope %.2f where the reasoning says %.2f. The curve there is a step: almost flat from 2^10 to 2^21 (%.2f to %.2f rounds) and rising after that, because nearly every prime needs 2 or 3 rounds. A second measure, chosen after seeing the step and therefore worth less: from the first class to the last the mean rises by %s rounds for k = %s, where the reasoning says %s. So the size of the effect is right within about a tenth in all four cases, and the straight line is a fair picture only when one round is short compared with the range of sizes.' % (s5['slope'], s5['predicted_slope'], s5['flat'][0], s5['flat'][1], ', '.join('%.2f' % d['rise'] for d in I['sources']), ', '.join(str(d['k']) for d in I['sources']), ', '.join('%.2f' % d['pred_rise'] for d in I['sources']))) if s5 else '') + '\n\n')
dp = I['deep']
s2 = dp['stage2']
r0 = s2['rows'][0]
DEEP = ('Primes of 7 rounds exist (computed; a search by construction, ' + L('inherit_deep.py') + ' → ' + L('inherit_deep_first.out') + ', ' + L(s2['file']) + '). Each neighbour of a prime of g + 1 rounds has a prime factor of g rounds or more (proved: otherwise that neighbour is lit completely by round g − 1, and N − 1 alone is already larger than √N). So candidates were built on the deep primes already known. From the %s primes r of 5 or 6 rounds below 4.2·10^9 (' % f(dp['deep5']) + L('inherit_deep_4e9.txt') + '), q = s · r ± 1 with even s up to %s gave %s primes of 6 rounds or more; from these, N = s · q ± 1 below %s gave %d primes of 7 rounds. The smallest is %s = %d × %s %s 1; the box of the page proves it in 7 rounds as well. This search finds examples. It does not show that this is the smallest prime of 7 rounds, and whether the number of rounds is unbounded is not proved.\n\n' % (f(dp['s1']), f(dp['six']), f(s2['bound']), s2['found'], f(r0[0]), r0[2], f(r0[3]), '+' if r0[4] == '+' else '−'))
s2p = I.get('second')
SECOND = (' A second C program, written separately from the description of the rule (' + L('inherit_second.c') + ' → ' + L('inherit_second_4e9.out') + '; primes and factors from one sieve, no Miller–Rabin test), gives the same number of primes for every round up to %s and the same list of %s primes of 5 and 6 rounds.' % (f(s2p['X']), f(s2p['deep']))) if s2p else ''
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

The runs, by three programs.

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

{UNSEEN}{SRC}{DEEP}What is known, and what we did not find. The two proof steps are Pocklington's (1914) and Morrison's (1975). Counting levels of such proofs is known in two forms: the classes of Erdős and Selfridge (all prime factors of p + 1, or of p − 1, down to 2 and 3; OEIS A005113 and A056637) and the Pratt tree (all prime factors of p − 1), whose height Ford, Konyagin and Luca (2010) bound for almost all primes and for which their model suggests about e · ln ln p. The count used here differs: both neighbours, and only as many helpers as it takes to pass √N. We did not find it in a search, and the local copy of the OEIS (31 October 2026) contains neither 2, 3, 19, 173, 21767, 1688173, 2223990533 nor the list of round 1. It is a small step from known things and may well be known.

What machine 23 is not. It is not a faster way to find primes: a round is heavy, because every dark wheel tries in every round, and a machine that may turn as much as it likes within one tick can let each wheel decide alone, as machine 20 does. For a large N the hard part is taking N − 1 and N + 1 apart. What the rounds measure is how deep a proof by these two steps has to go.

On the page ({L('page_check_20_23.txt')}; the page's own code in headless Chrome with 5,000 wheels): machine 20 has 669 lamps at tick 12, machine 21 at tick 19, machine 23 after round 3, machine 22 at tick 5,000 with {f(wk['count_primes'])} primes counted; no mistake on any of them. The page code of the four machines was compared with Python ({LG('early_node.js')}, {L('early_js_check.py')} → {L('early_js_check.txt')}: {E['js']} comparisons, 0 differ; {LG('inherit_node.js')}, {L('inherit_js_check.py')} → {L('inherit_js_check.txt')}: rounds, helpers, witnesses and the division for {f(I['js'])} odd numbers, 0 differ).

'''
before('**What this is and is not.**', text)
rep("the chain of wheels and the test with a receipt are literature and their runs are computed)", "the chain of wheels and the test with a receipt are literature and their runs are computed; the four machines that find primes before their tick are literature with computed runs, and the rule for the rounds of machine 23 is a heuristic)")
open('section20.md', 'w', encoding='utf-8').write(s)
print(len(s))
