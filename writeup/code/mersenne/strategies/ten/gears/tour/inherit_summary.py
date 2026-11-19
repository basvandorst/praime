"""Collects the results of machine 23 (wheels that light each other) into inherit_sum.json and inherit_sum.txt.
Every number is read from the output of the program that computed it."""
import json, re, math, os
from sympy import isprime, primerange
from early import lucy
from inherit import run

out = {}
lines = []


def say(s):
    lines.append(s)
    print(s)


P = json.load(open('inherit.json'))
both = P['n - 1 and n + 1']
minus = P['n - 1 only']
assert both['ok'] and minus['ok'] and P['composites'][1] == 0
out['python'] = P
say('proofs run for every odd number up to %d (Python): the lit numbers are exactly the primes, with both neighbours and with n - 1 alone; %d composites with full helpers: %d passed' % (both['top'], P['composites'][0], P['composites'][1]))


def read_c(name):
    t = open(name).read()
    X = int(re.search(r'X (\d+) side', t).group(1))
    cnt = [(int(a), int(b), int(c)) for a, b, c in re.findall(r'depth (\d+): (\d+) primes, first (\d+)', t)]
    mean = [(int(a), int(b), float(c)) for a, b, c in re.findall(r'primes in \[2\^(\d+), 2\^\d+\): (\d+), mean depth ([\d.]+)', t)]
    mean = [m for m in mean if 2 ** (m[0] + 1) - 1 <= X]  # the last class of a run is cut off by X: leave it out
    assert 'done' in t
    return X, cnt, mean


X, cnt, mean = read_c('inherit_c_4e9.out')
V, C, moves, pr = lucy(X)
assert sum(c[1] for c in cnt) == C[0]
for d, n, first in cnt:
    if d < len(both['first']):
        assert both['first'][d][1] == first
py = dict(both['cnt'])
small = read_c  # the C program on the range of the Python run is compared in final.sh (inherit 200000)
out['c'] = dict(X=X, cnt=cnt, mean=mean, primes=C[0])
say('C program up to %d: %d primes (equal to the count of the counting wheels); by rounds: %s' % (X, C[0], ', '.join('%d: %d' % (d, n) for d, n, f in cnt)))
say('  smallest prime that needs 1, 2, ... rounds: %s; most rounds: %d' % (', '.join(str(f) for d, n, f in cnt[1:]), cnt[-1][0]))
Xm, cntm, meanm = read_c('inherit_c_minus_1e9.out')
out['c_minus'] = dict(X=Xm, cnt=cntm, mean=meanm)
say('with n - 1 alone, up to %d: most rounds %d; smallest prime of each round: %s' % (Xm, cntm[-1][0], ', '.join(str(f) for d, n, f in cntm[1:])))
say('  mean rounds for primes between 2^28 and 2^29: %.2f with both neighbours, %.2f with n - 1 alone' % (dict((a, c) for a, b, c in mean)[28], dict((a, c) for a, b, c in meanm)[28]))
out['fit_classes'] = [mean[0][0] if mean[0][0] >= 10 else 10, mean[-1][0]]
big = [json.loads(l) for l in open('inherit_big.out') if l.startswith('{')]
out['big'] = big
for b in big:
    say('  random primes of %d digits: %d tried, mean %.2f rounds, most %d; by rounds %s' % (b['digits'], b['count'], b['mean'], b['worst'][0], b['hist']))
# least squares: mean rounds against ln ln N, classes from 2^10 on and the samples
pts = [(math.log(math.log(1.5 * 2 ** a)), c) for a, b, c in mean if a >= 10] + [(math.log(math.log(5.5 * 10 ** (b['digits'] - 1))), b['mean']) for b in big if b['digits'] >= 12]
n = len(pts)
sx = sum(p[0] for p in pts); sy = sum(p[1] for p in pts)
sxx = sum(p[0] ** 2 for p in pts); sxy = sum(p[0] * p[1] for p in pts)
slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
icpt = (sy - slope * sx) / n
resid = max(abs(p[1] - (icpt + slope * p[0])) for p in pts)
out['fit'] = dict(slope=round(slope, 3), intercept=round(icpt, 3), points=n, max_resid=round(resid, 3))
say('  straight line through %d points (mean rounds against ln ln N): rounds = %.2f + %.2f ln ln N, largest distance %.2f' % (n, icpt, slope, resid))
# the same points with the slope fixed at 1 / ln 4 (both neighbours) and 1 / ln 2 (n - 1 alone): how far the constant moves
c4 = [p[1] - p[0] / math.log(4) for p in pts]
pm = [(math.log(math.log(1.5 * 2 ** a)), c) for a, b, c in meanm if 10 <= a <= 28]
c2 = [p[1] - p[0] / math.log(2) for p in pm]
out['law'] = dict(both=[round(min(c4), 3), round(max(c4), 3), round(sum(c4) / len(c4), 3)], minus=[round(min(c2), 3), round(max(c2), 3), round(sum(c2) / len(c2), 3)], minus_last_slope=round((pm[-1][1] - pm[-4][1]) / (pm[-1][0] - pm[-4][0]), 2))
say('  with the slope fixed at 1/ln 4 = %.4f: rounds - log4(ln N) lies between %.2f and %.2f for all %d points (mean %.2f)' % (1 / math.log(4), min(c4), max(c4), n, sum(c4) / len(c4)))
say('  n - 1 alone, classes 2^10 to 2^28, slope fixed at 1/ln 2: rounds - log2(ln N) lies between %.2f and %.2f; slope over the last three classes %.2f' % (min(c2), max(c2), out['law']['minus_last_slope']))
B = open('inherit_big_check.txt').read()
m = re.search(r'(\d+) primes up to (\d+): (\d+) differ', B)
assert m.group(3) == '0' and '2223990533: 6 rounds' in B
out['topdown'] = int(m.group(1))
say('the rule from the top down (PARI) against the run from the bottom up: %s primes, 0 differ; 2,223,990,533 needs 6 rounds there too' % m.group(1))
H = open('inherit_hist.txt').read()
m = re.search(r'equal: True \((\d+) primes\)', H)
m2 = re.search(r'every prime from 3 to (\d+)', H)
m3 = re.search(r'the (\d+) primes of 5 and 6 rounds .*?: ([\d, ]+)\): rounds from the top down differ for (\d+)', H)
assert m and m3.group(3) == '0'
out['hist'] = dict(top=int(m2.group(1)), primes=int(m.group(1)), six=[int(x) for x in m3.group(2).split(', ')], listed=int(m3.group(1)))
say('third count (PARI, top down) for every prime up to %s: histogram equal to the C program (%s primes); the primes of 6 rounds up to 4.2e9 are %s' % (m2.group(1), m.group(1), m3.group(2)))
S2 = None
if os.path.exists('inherit_second_4e9.out') and 'total primes' in open('inherit_second_4e9.out').read():
    t2 = open('inherit_second_4e9.out').read()
    c2 = [(int(a), int(b), int(c)) for a, b, c in re.findall(r'round (\d+): (\d+) primes, first (\d+)', t2)]
    deep2 = sorted(int(b) for a, b in re.findall(r'deep (\d+) (\d+)', t2))
    deep1 = sorted(int(l.split()[1]) for l in open('inherit_deep_4e9.txt'))
    assert c2 == [tuple(x) for x in cnt] and deep2 == deep1, 'second program differs'
    S2 = dict(X=int(re.search(r'X2 (\d+)', t2).group(1)), deep=len(deep2))
    say('second C program, written separately (primes from a sieve): up to %d the primes per round and the %d primes of 5 or 6 rounds are the same' % (S2['X'], S2['deep']))
out['second'] = S2
S8 = None
if os.path.exists('inherit_second_8e9.out') and 'total primes' in open('inherit_second_8e9.out').read():
    t8 = open('inherit_second_8e9.out').read()
    c8 = [(int(a), int(b), int(c)) for a, b, c in re.findall(r'round (\d+): (\d+) primes, first (\d+)', t8)]
    X8 = int(re.search(r'X2 (\d+)', t8).group(1))
    V8, C8, mv8, pr8 = lucy(X8)
    assert sum(c[1] for c in c8) == C8[0] == int(re.search(r'total primes (\d+)', t8).group(1))
    six8 = sorted(int(b) for a, b in re.findall(r'deep (\d+) (\d+)', t8) if int(a) >= 6)
    LM = open('inherit_lemma.txt').read()
    lem = [int(x) for x in re.findall(r'(\d+) \(\d+ rounds\)', LM)]
    assert lem == six8 and c8[-1][0] == 6, (lem, six8)
    S8 = dict(X=X8, cnt=c8, primes=C8[0], six=six8, lemma_cands=int(re.search(r': (\d+), of which', LM).group(1)))
    d8 = sum(1 for l in open('inherit_deep_8e9.txt'))
    assert d8 == sum(c[1] for c in c8 if c[0] >= 5)

    def rd(name):
        t = open(name).read()
        a = re.search(r'from the (\d+) primes of 5 or more rounds up to (\d+): (\d+) candidates s r \+- 1 up to (\d+), primes of 6 or more rounds: (\d+) \(most rounds 6\)', t)
        b = re.search(r'from these: (\d+) candidates s q \+- 1 up to (\d+), primes of 7 or more rounds: none', t)
        assert a and b, name
        return dict(list=int(a.group(1)), B=int(a.group(2)), cands=int(a.group(3)), top6=int(a.group(4)), six=int(a.group(5)), cands7=int(b.group(1)), top7=int(b.group(2)))

    ch4, ch8 = rd('inherit_lemma_chain_4e9.txt'), rd('inherit_lemma_chain_8e9.txt')
    assert ch4['list'] == S2['deep'] and ch4['B'] == S2['X'] and ch8['list'] == d8 and ch8['B'] == X8
    S8.update(deep=d8, chain4=ch4, chain8=ch8, no7=ch8['top7'], no7_both=ch4['top7'])
    say('  the lemma used twice: from the list up to %d: %d primes of 6 rounds up to %d and no prime of 7 rounds up to %d; from the list up to %d: %d primes of 6 rounds up to %d and no prime of 7 rounds up to %d' % (ch4['B'], ch4['six'], ch4['top6'], ch4['top7'], ch8['B'], ch8['six'], ch8['top6'], ch8['top7']))
if S8 and os.path.exists('inherit_third_17e9.out') and 'total primes' in open('inherit_third_17e9.out').read():
    t3 = open('inherit_third_17e9.out').read()
    c3 = [(int(a), int(b), int(c)) for a, b, c in re.findall(r'round (\d+): (\d+) primes, first (\d+)', t3)]
    X3 = int(re.search(r'X2 (\d+)', t3).group(1))
    V3, C3, mv3, pr3 = lucy(X3)
    assert sum(c[1] for c in c3) == C3[0] and c3[-1][0] == 6
    six3 = sorted(int(b) for a, b in re.findall(r'deep (\d+) (\d+)', t3) if int(a) >= 6)
    lem11 = [int(x) for x in re.findall(r'(\d+) \(6 rounds\)', open('inherit_lemma_2.txt').read())]
    d3 = [int(l.split()[1]) for l in open('inherit_deep_17e9.txt')]
    d8l = [int(l.split()[1]) for l in open('inherit_deep_8e9.txt')]
    assert six3 == lem11 and len(six3) == S8['chain4']['six'] and len(d3) == sum(c[1] for c in c3 if c[0] >= 5) and d3[:len(d8l)] == d8l
    ch3 = rd('inherit_lemma_chain_17e9.txt')
    assert ch3['list'] == len(d3) and ch3['B'] == X3 and ch3['six'] >= S8['chain8']['six']
    cl3 = [(int(a), int(b), float(c)) for a, b, c in re.findall(r'class 2\^(\d+): (\d+) primes, mean ([\d.]+)', t3)]
    for a, b, c in cl3:
        if a <= mean[-1][0] and a >= 10:
            assert abs(dict((x, z) for x, y, z in mean)[a] - c) < 6e-5, a
    drift = [(a, round(c - math.log(math.log(1.5 * 2 ** a)) / math.log(4), 3)) for a, b, c in cl3 if 29 <= a and 2 ** (a + 1) - 1 <= X3]
    newc = [dict(a=a, primes=b, mean=c, pred=round(math.log(math.log(1.5 * 2 ** a)) / math.log(4) + round(out['law']['both'][2], 2), 3)) for a, b, c in cl3 if a > mean[-1][0] and 2 ** (a + 1) - 1 <= X3]
    S8.update(third=dict(X=X3, cnt=c3, primes=C3[0], six=six3, deep=len(d3), chain=ch3, classes=newc, drift=drift), no7=ch3['top7'])
    for e in newc:
        say('  class not used for the rule: every prime from 2^%d to 2^%d (%d primes): mean %.3f rounds, the rule says %.3f' % (e['a'], e['a'] + 1, e['primes'], e['mean'], e['pred']))
    say('third program (the second with a smaller table) up to %d: %d primes (equal to the counting recursion), by rounds %s; the %d primes of 6 rounds are those of the lemma; with its %d primes of 5 or more rounds: %d primes of 6 rounds up to %d, no prime of 7 rounds up to %d' % (X3, C3[0], ', '.join('%d: %d' % (d, n) for d, n, f in c3), len(six3), len(d3), ch3['six'], ch3['top6'], ch3['top7']))
out['second8'] = S8
X = open('inherit_box_check.txt').read()
m = re.search(r'box of machine 23: (\d+) primes .*?failed: (\d+); primes where the page could not take a neighbour apart and says so: (\d+); not finished by the page: (\d+)', X)
assert m.group(2) == '0' and m.group(4) == '0'
mc = re.search(r'helpers of 2223990533 on the page, largest first: ([\d, ]+)', X)
md = re.search(r'they have (\d+) to (\d+) digits', X)
out['boxcheck'] = dict(primes=int(m.group(1)), partial=int(m.group(3)), chain=[int(x) for x in mc.group(1).split(', ')], digits=[int(md.group(1)), int(md.group(2))])
say('box of machine 23 against PARI: %s primes, 0 differ, %s where the page says it could not take a neighbour apart' % (m.group(1), m.group(3)))
un = []
by = {}
for name in ('inherit_big_65.out', 'inherit_big_70.out', 'inherit_big_75.out', 'inherit_big_75b.out', 'inherit_big_80.out', 'inherit_big_85.out'):
    if os.path.exists(name) and os.path.getsize(name):
        b = json.loads(open(name).read().strip().split('\n')[-1])
        e = by.setdefault(b['digits'], dict(digits=b['digits'], files=[], hist={}, most=0))
        e['files'].append(name)
        e.setdefault('parts', []).append((name, b['count'], dict(b['hist'])))
        for k_, v_ in b['hist']:
            e['hist'][k_] = e['hist'].get(k_, 0) + v_
        e['most'] = max(e['most'], b['worst'][0])
for dg in sorted(by):
    e = by[dg]
    nn = sum(e['hist'].values())
    mu = sum(k_ * v_ for k_, v_ in e['hist'].items()) / nn
    se = (sum(v_ * (k_ - mu) ** 2 for k_, v_ in e['hist'].items()) / (nn - 1) / nn) ** 0.5
    pred = math.log(math.log(5.5 * 10 ** (dg - 1))) / math.log(4) + round(out['law']['both'][2], 2)
    parts = []
    for name, n1, h1 in e['parts']:
        m1 = sum(int(k_) * v_ for k_, v_ in h1.items()) / n1
        s1 = (sum(v_ * (int(k_) - m1) ** 2 for k_, v_ in h1.items()) / (n1 - 1) / n1) ** 0.5
        parts.append(dict(file=name, count=n1, mean=round(m1, 3), se=round(s1, 3), z=round((m1 - pred) / s1, 1)))
    un.append(dict(digits=dg, count=nn, mean=round(mu, 3), se=round(se, 3), pred=round(pred, 3), most=e['most'], hist=sorted(e['hist'].items()), files=e['files'], z=round((mu - pred) / se, 1), parts=parts))
    say('  size not used for the rule: %d random primes of %d digits, mean %.3f rounds (standard error %.3f), the rule says %.3f: %.1f standard errors; most %d' % (nn, dg, mu, se, pred, (mu - pred) / se, e['most']))
out['unseen'] = un
src = []
for k in (1, 2, 3, 5):
    name = 'inherit_sources_%d.out' % k
    if os.path.exists(name) and os.path.getsize(name):
        d = json.loads(open(name).read().strip())
        d['rise'] = round(d['rows'][-1][3] - d['rows'][0][3], 3)
        d['pred_rise'] = round((d['rows'][-1][1] - d['rows'][0][1]) / math.log(2 * k), 3)
        d['flat'] = [d['rows'][0][3], d['rows'][5][3]]
        src.append(d)
        say('  %d source%s of helpers: slope of mean rounds against ln ln N %.2f, the reason says 1/ln %d = %.2f (%d size classes; rounds - log_%d(ln N) between %.2f and %.2f)' % (k, '' if k == 1 else 's', d['slope'], 2 * k, d['predicted_slope'], len(d['rows']), 2 * k, d['const_min'], d['const_max']))
        say('     rise of the mean from the first class to the last: %.2f rounds, the reason says %.2f' % (d['rise'], d['pred_rise']))
out['sources'] = src
two = [d for d in src if d['k'] == 2][0]
cl = {a: c for a, b, c in mean}
for name, x, n_, m_ in two['rows']:
    if name.startswith('2^'):
        assert abs(cl[int(name[2:4])] - m_) < 6e-5, (name, m_)
D1 = open('inherit_deep_first.out').read()
m = re.search(r'primes of 5 or more rounds below 4.2e9: (\d+)', D1)
m1 = re.search(r'stage 1: (\d+) candidates s r \+- 1 with s <= (\d+), (\d+) primes, (\d+) of them need 6 rounds or more', D1)
best = None
for name in ('inherit_deep2.out', 'inherit_deep_wide2.out'):
    if os.path.exists(name) and 'done' in open(name).read():
        D2 = open(name).read()
        m2 = re.search(r'stage 2 below (\d+): (\d+) primes q of 6 rounds or more, (\d+) candidates s q \+- 1, (\d+) primes, (\d+) need 7 rounds or more', D2)
        rows = [(int(a), int(b), int(c), int(d), e) for a, b, c, d, e in re.findall(r'  (\d+) needs (\d+) rounds \((\d+) \* (\d+) ([+-]) 1\)', D2)]
        cand = dict(file=name, bound=int(m2.group(1)), q=int(m2.group(2)), cands=int(m2.group(3)), primes=int(m2.group(4)), found=len(set(r[0] for r in rows)) if int(m2.group(5)) <= 12 else int(m2.group(5)), rows=rows)
        if best is None or rows[0][0] < best['rows'][0][0]:
            best = cand
assert best and best['rows'][0][1] == 7
out['deep'] = dict(deep5=int(m.group(1)), s1=int(m1.group(2)), cands1=int(m1.group(1)), primes1=int(m1.group(3)), six=int(m1.group(4)), stage2=best)
r0 = best['rows'][0]
if os.path.exists('inherit_deep_wide2.out') and 'done' in open('inherit_deep_wide2.out').read():
    W1 = open('inherit_deep_wide.out').read()
    W2 = open('inherit_deep_wide2.out').read()
    w1 = re.search(r'stage 1: (\d+) candidates s r \+- 1 with s <= \d+, (\d+) primes, (\d+) of them need 6 rounds or more', W1)
    w2 = re.search(r'stage 2 below (\d+): (\d+) primes q of 6 rounds or more, (\d+) candidates s q \+- 1, (\d+) primes', W2)
    wmin = min(int(a) for a in re.findall(r'  (\d+) needs 7 rounds', W2))
    assert int(w1.group(3)) == int(w2.group(2)) and wmin == r0[0]
    out['deep']['wide'] = dict(cands=int(w1.group(1)), primes=int(w1.group(2)), six=int(w1.group(3)), bound=int(w2.group(1)), cands2=int(w2.group(3)), smallest=wmin)
    say('  wider first stage: %s candidates, %s primes of 6 rounds or more; below %s the smallest prime of 7 rounds found is the same' % (w1.group(1), w1.group(3), w2.group(1)))
say('search by construction: %s primes of 5 or 6 rounds below 4.2e9; q = s r +- 1 with s <= %s: %s primes of 6 rounds or more; N = s q +- 1 below %d: %d primes of 7 rounds; the smallest found is %d = %d * %d %s 1 (not shown to be the smallest of all)' % (m.group(1), m1.group(2), m1.group(4), best['bound'], best['found'], r0[0], r0[2], r0[3], r0[4]))
X7 = open('inherit_box_check.txt').read()
assert 'four primes of 7 rounds' in X7 and out['boxcheck']['chain'][0] == 555997633
J = open('inherit_js_check.txt').read()
m = re.search(r'(\d+) odd numbers, .*?, (\d+) differ', J)
assert m.group(2) == '0'
out['js'] = int(m.group(1))
say('page code against Python: %s odd numbers, 0 differ' % m.group(1))
d = run(5000, 2)
MX = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127}
Pp = list(primerange(2, 251))
TW = [p for p in Pp if isprime(p - 2) or isprime(p + 2)]
out['row250'] = [len(Pp), max(d[p] for p in Pp), len([p for p in Pp if p in MX]), max(d[p] for p in Pp if p in MX), len(TW), max(d[p] for p in TW)]
c5 = {}
for p in d:
    c5[d[p]] = c5.get(d[p], 0) + 1
out['upto5000'] = sorted(c5.items())
r1 = sorted(p for p in d if d[p] == 1)
proth = [p for p in r1 if ((p - 1) & -(p - 1)) ** 2 > p]
out['round1'] = dict(count=len(r1), proth=len(proth), first=r1[:14])
assert all((max((p - 1) & -(p - 1), (p + 1) & -(p + 1))) ** 2 > p for p in r1)
say('250 wheels: all %d primes lit from round %d, the %d Mersenne exponents from round %d, the %d twin primes from round %d; up to 5,000 by rounds: %s; round 1 up to 5,000: %d primes, %d of them of the form k 2^m + 1 with k < 2^m' % (tuple(out['row250']) + (out['upto5000'], len(r1), len(proth))))
json.dump(out, open('inherit_sum.json', 'w'))
open('inherit_sum.txt', 'w').write('\n'.join(lines) + '\n')
