"""The three-cell wheel (fold.c -three): a cubic cousin of the two-cell wheel, built by folding Pascal's ring.
Wheel sizes r: primes r = 1 (mod 3) below 256 for which 2 is not a cube mod r. For n take the first r that does not divide n
and for which n is not a cube mod r; the wheel has three cells (the sums of x^k over the three classes of the cubes mod r)
and one turn multiplies by u = product of (x^h + 1) over the cubes h. Home: after n turns the cells show the start, moved round by the class of n.
Checks: every prime comes home (proved: Frobenius), the C program against plain Python on the full ring, and the runs on numbers that are not prime."""
import subprocess, random, re, sys, json
from sympy import primerange, isprime
from fold_check import fold_vec, classes
def run(nums):
    o = subprocess.run(['./fold', '-three'], input='\n'.join(map(str, nums)) + '\n', capture_output=True, text=True)
    return set(int(l.split()[4]) for l in o.stdout.split('\n') if l.startswith('three')), o.stderr
sizes = [int(x) for x in re.search(r'wheel sizes:([ \d]+)', run([9])[1]).group(1).split()]
def wheel(n):
    for r in sizes:
        if n % r == 0:
            if n == r: continue
            return 0
        if classes(r, 3)[n % r] != 0: return r
    return None
random.seed(5)
ps = list(primerange(3, 300000)) + [18446744073709551557, 18446744073709551533, 18446744073709551521]
home, err = run(ps)
print('wheel sizes:', sizes)
print('primes from 3 to 300000 and three primes just below 2^64:', len(ps), '| came home:', len(home), '| did not:', [p for p in ps if p not in home][:10])
big = [int(x) for x in open('sample64.txt')]; psp = [int(x) for x in open('psp2_to1e10.txt')]
pick = random.sample(big, 150) + random.sample(psp, 150) + random.sample(ps, 100) + [n for n in range(9, 3000, 2)]
h2, _ = run(pick); cmpd = wrong = 0
for n in pick:
    r = wheel(n)
    ref = bool(r) and fold_vec(n, r, 3)[0]
    cmpd += 1
    if ref != (n in h2): wrong += 1
print('compared with plain Python on the full ring:', cmpd, 'numbers | different verdicts:', wrong)
R = {'sizes': sizes}
for tag in sys.argv[1:]:
    e = open('three_%s.err' % tag).read(); o = open('three_%s.out' % tag).read().strip()
    m = re.search(r'numbers read (\d+), last (\d+), came home (\d+), turned away \(divisible by a wheel size\) (\d+), no wheel (\d+)', e)
    R[tag] = dict(zip(('tested', 'last', 'home', 'away', 'nowheel'), map(int, m.groups()))); R[tag]['use'] = {a: int(b) for a, b in re.findall(r'ring of (\d+) cells used (\d+) times', e)}
    assert (not o) == (R[tag]['home'] == 0)
    print('%s: %s numbers, %d came home, %s turned away (divisible by a wheel size), %d with no wheel; wheels used: %s' % (tag, format(R[tag]['tested'], ','), R[tag]['home'], format(R[tag]['away'], ','), R[tag]['nowheel'], ', '.join('%s: %s' % (a, format(b, ',')) for a, b in R[tag]['use'].items())))
json.dump(R, open('three_check.json', 'w'), indent=1)
