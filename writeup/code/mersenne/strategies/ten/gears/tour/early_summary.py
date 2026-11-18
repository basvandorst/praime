"""Collects the results of machines 20, 21 and 22 into early_sum.json and early_sum.txt.
Every number is read from the output of the program that computed it, and the counts of the squaring wheels
per power of two (C program norm2, early_dyadic.sh) are compared with the counting wheels (early.lucy)."""
import json, re
from sympy import isprime, primerange, primepi
from early import lucy

out = {}
lines = []


def say(s):
    lines.append(s)
    print(s)


E = json.load(open('early.json'))
assert E['ladder']['bad'] == 0 and E['ladder']['badt'] == 0 and E['self']['wrong'] == 0
out['early'] = E
# squaring wheels: numbers that came home in [2^k, 2^(k+1)), k = 1..31
dy = {}
for m in re.finditer(r'k (\d+) norm-2 wheel, odd n in \[(\d+), (\d+)\): came home (\d+)', open('early_dyadic.out').read()):
    dy[int(m.group(1))] = int(m.group(4))
assert sorted(dy) == list(range(1, 32))
cum = 1
table = []
bad = 0
for T in range(1, 32):
    cum += dy[T]
    V, C, moves, pr = lucy(2 ** (T + 1) - 1)
    if C[0] != cum:
        bad += 1
    table.append((T, 2 ** (T + 1), cum))
assert bad == 0
out['dyadic'] = table
say('squaring wheels: lamps lit by tick T = numbers below 2^(T+1) that came home (C program, every odd number), T = 1..31; equal to the count of the counting wheels for every T: differences %d' % bad)
say('  tick 10: %d, tick 12: %d, tick 20: %d, tick 30: %d, tick 31: %d' % (table[9][2], table[11][2], table[19][2], table[29][2], table[30][2]))
out['pi_2_64'] = 425656284035217743
say('  tick 63 (every wheel below 2^64): 425,656,284,035,217,743 primes (literature: OEIS A007053), given the list of machine 19')
# rows of the table "most primes at home" for 250 wheels
MX = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127}
P = list(primerange(2, 251))
TW = [p for p in P if isprime(p - 2) or isprime(p + 2)]
gear_tick = {}
for r in E['self']['rows']:
    gear_tick[r[0]] = r[3]


def sf_tick(k):
    return next(t for t in sorted(gear_tick) if gear_tick[t] > k)


for key, tick in (('sq', lambda k: k.bit_length() - 1), ('sf', sf_tick), ('cnt', lambda k: k)):
    out['row250_' + key] = [len(P), max(tick(p) for p in P), len([p for p in P if p in MX]), max(tick(p) for p in P if p in MX), len(TW), max(tick(p) for p in TW)]
    say('250 wheels, machine %s: all %d primes at home from tick %d, all %d Mersenne exponents from tick %d, all %d twin primes from tick %d' % ((key,) + tuple(out['row250_' + key])))
# work per prime with N = 5000 wheels (x = N^2 for the counting wheels)
N = 5000
pr = list(primerange(2, N + 1))
ladder_ticks = sum(k.bit_length() - 1 for k in range(3, N + 1, 2))
strikes = sum((N - p * p) // p + 1 for p in pr if p * p <= N)
V, C, moves, _ = lucy(N * N)
out['work'] = dict(N=N, primes=len(pr), ladder_ticks=ladder_ticks, strikes=strikes, looked=N - 1, count_moves=moves, count_primes=C[0])
say('work for 5,000 wheels: squaring wheels %d ticks of an odd wheel (%.1f per prime found); gears from lamps %d strikes and %d numbers looked at (%.1f per prime); counting wheels %d moves for %d primes counted (%.3f per prime)' % (ladder_ticks, ladder_ticks / len(pr), strikes, N - 1, (strikes + N - 1) / len(pr), moves, C[0], moves / C[0]))
J = open('early_js_check.txt').read()
m = re.search(r'(\d+) comparisons, (\d+) differ', J)
assert m.group(2) == '0' and 'differ from the true number of primes or from the moves of the Python twin: 0' in J
out['js'] = int(m.group(1))
m = re.search(r'largest x: pi\((\d+)\) = (\d+) in (\d+) moves', J)
out['box'] = [int(x) for x in m.groups()]
say('page code against Python: %d comparisons, 0 differ; the box of machine 22 counted %d primes up to %d in %d moves' % (out['js'], out['box'][1], out['box'][0], out['box'][2]))
json.dump(out, open('early_sum.json', 'w'))
open('early_sum.txt', 'w').write('\n'.join(lines) + '\n')
