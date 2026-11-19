"""Second, separate check of mexp_gaps.json: recomputes the gaps, the size statistics and the counts of gaps that are
multiples of m with sympy and plain integers, and compares model A with its exact value.
Model A draws each exponent as a random prime of about its size. For m not dividing any exponent >= 61, the class of a
random large prime mod m is (close to) uniform over the phi(m) classes prime to m, so two independent ones agree with
probability 1/phi(m): the expected count is (number of gaps) / phi(m)."""
import json, math
from sympy import isprime, totient, primerange

J = json.load(open('mexp_gaps.json'))
E = J['exponents']
bad = 0
assert len(E) == 52 and all(isprime(p) for p in E)
# every listed exponent really gives a Mersenne prime is literature; here: 2^p - 1 has no factor 2kp + 1 for small k
for p in E[:20]:
    if p > 2:
        for k in range(1, 200):
            q = 2 * k * p + 1
            if q < 2 ** p - 1 and isprime(q) and pow(2, p, q) == 1:
                bad += 1
gaps = [b - a for a, b in zip(E, E[1:])]
bad += gaps != J['gaps']
ex = E[:51]
for start in (5, 127):
    seq = [p for p in ex if p >= start]
    r = [math.log(b / a, 2) for a, b in zip(seq, seq[1:])]
    s = J['sizes_from_%d' % start]
    bad += abs(sum(r) / len(r) - s['mean']) > 1e-9 or len(r) != s['n']
    bad += abs(s['model'] - math.exp(-0.5772156649015329)) > 1e-12
seq = [p for p in ex if p >= 61]
pairs = len(seq) - 1
bad += pairs != J['mods']['pairs']
lines = []
for row in J['mods']['rows']:
    m = row['m']
    real = sum(1 for a, b in zip(seq, seq[1:]) if (b - a) % m == 0)
    exact = pairs / int(totient(m))
    bad += real != row['real']
    ok = abs(row['A'] - exact) < 0.6
    bad += not ok
    lines.append('m = %d: %d gaps are multiples of m; model A (simulated) %.2f, exactly %d/phi(%d) = %.2f' % (m, real, row['A'], pairs, m, exact))
c4 = sum(1 for p in seq if p % 4 == 1)
bad += c4 != J['class4']['real']
cnt = [sum(1 for p in ex if 2 ** k <= p < 2 ** (k + 1)) for k in range(2, 26)]
bad += cnt != J['doublings']['counts']
for w in J['wheels']['rows']:
    bad += w['real'] != sum(1 for a, b in zip(seq, seq[1:]) if (b - a) % w['m'] == 0)
print('\n'.join(lines))
print('gaps, sizes, counts per doubling, %d remainder rows and %d wheel counts recomputed: %d differ' % (len(J['mods']['rows']), len(J['wheels']['rows']), bad))
