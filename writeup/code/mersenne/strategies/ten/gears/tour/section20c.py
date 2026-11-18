"""Third round of section 20: reads section20.md as written by section20.py and section20b.py and adds the chain of wheels (machine 17).
Every number comes from chain.json, which chain_summary.py builds from the output files of the test scripts."""
import json
CH = json.load(open('chain.json'))
f = lambda n: format(n, ',')
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


modes = [('all sizes of wheel, 97 curve families', 'n − 1, n + 1 and 97 curve families'), ('curves only, 97 families', '97 curve families only'), ('all sizes of wheel, 9 curve families', 'n − 1, n + 1 and 9 curve families'), ('n - 1 and n + 1 only', 'n − 1 and n + 1 only')]
digs = sorted(set(r['digits'] for r in CH['big']))


def cell(d, m):
    r = [x for x in CH['big'] if x['digits'] == d and x['mode'].strip() == m]
    if not r:
        return 'not run'
    r = r[0]
    t = '%d of %d' % (r['proved'], r['primes'])
    if r['timeout']:
        t += ' (%d stopped after %d s)' % (r['timeout'], r['tlim'])
    return t


def med(d):
    r = [x for x in CH['big'] if x['digits'] == d and x['mode'].strip() == modes[0][0]][0]
    return r['steps'], r['median_s']


rows = '\n'.join('| %d | %s | %s | %s |' % (d, ' | '.join(cell(d, m) for m, _ in modes), med(d)[0], med(d)[1]) for d in digs)
sm = CH['small']
adv = CH['adv']
npq = sum(sum(adv[k].values()) for k in adv if k.startswith('pq'))
j = {r['digits']: r for r in CH['js'] if 'digits' in r}
jr = [r for r in CH['js'] if 'digits' not in r]
text = f'''**A chain of wheels (the three theorems are literature; the runs are computed).** Machines 14, 15 and 16 are one idea with three sources of wheel sizes, and together they give a machine that is proved, fast, and has worked for every number tried. Let n be the number and f, q primes.

1. One cell (Pocklington). Suppose a^(n−1) ≡ 1 (mod n) and gcd(a^((n−1)/f) − 1, n) = 1 for every prime f of a set, and let F be the product of the full powers of these f in n − 1. For a prime p dividing n, the order of a mod p divides n − 1 and not (n − 1)/f, so it contains the full power of f; it divides p − 1. Hence p ≡ 1 (mod F). If F² > n, every prime factor of n exceeds √n, and n is prime.
2. Two cells (the n + 1 partner; Morrison). Let c be coprime to n, u = a + √c with a² − c coprime to n, and β = u²/(a² − c), which is u divided by its conjugate and has norm 1. Suppose β^(n+1) = 1 and β^((n+1)/f) − 1 has a norm coprime to n for every prime f of a set, with F the product of their full powers in n + 1. Mod a prime p dividing n the two cells are either two copies of the remainders mod p (c is a square mod p; orders divide p − 1) or the field of p² elements, where the elements of norm 1 form a group of order p + 1. As before the order of β contains F, so F divides p − 1 or p + 1 and p ≥ F − 1. If (F − 1)² > n, n is prime.
3. A curve (Goldwasser and Kilian). Let gcd(n, 6) = 1, let E be y² = x³ + ax + b with gcd(4a³ + 27b², n) = 1, P a point of E mod n, and m = c · q with q > (n^(1/4) + 1)². Add points by the chord-and-tangent rule, where every division is by a number coprime to n. If m · P is the point at infinity and c · P is not, n is prime. Proof: mod a prime p dividing n every step of the same computation is a valid step on E over the remainders mod p, so the order of P there divides m and not c; q is prime, so q divides that order. The curve has at most (√p + 1)² points (Hasse), so (√p + 1)² ≥ q > (n^(1/4) + 1)² and p > √n.

In each case n is prime if the large prime q (or f) is, and that number is smaller: a chain, ended by trial division below 10^10.

Where the curve sizes come from (Atkin and Morain). Machine 16 counts the points of a curve, which takes about n steps. For a negative discriminant D with 4n = t² + |D| s² the curves with complex multiplication by D have n + 1 − t and n + 1 + t points (for D = −4 also n + 1 ± 2s, for D = −3 also n + 1 ± (t ± 3s)/2); t and s follow from a square root of D mod n by Cornacchia's algorithm. The j-invariants of these curves are the roots mod n of the class polynomial of D. We use the {CH['families']} fundamental discriminants with class number at most 4 ({L('classpoly_gen.py')} with PARI → {L('classpoly.json')}; {L('classpoly_check.py')} recomputes every polynomial from the j-function with mpmath: {CH['families']} checked, 0 differ). From a root j the curve is y² = x³ − 3kx + 2k with k = j/(j − 1728); for a random x0 and z = x0³ − 3k x0 + 2k the curve y² = x³ − 3k z² x + 2k z³ contains the point (z x0, z²), so no square root is needed, and z decides which of the two sizes it has. For D = −3 and −4 the curve is y² = x³ + b or y² = x³ + ax through a random point.

The search ({L('chain_wheel.py')}; the page's version is the same in JavaScript). (i) n must pass the doubling wheel (the strong test to base 2). (ii) List every size on offer: n − 1, n + 1 and the curve sizes. Divide out the primes below 10^5; keep a size when what remains is a probable prime q that is large enough. (iii) Take the candidates with the smallest q first; prove q by the same procedure, remembering the numbers that failed; then find a, c or the curve and the point, and run the check of the theorem. (iv) If nothing works, split the remaining parts with Pollard's rho (20,000 steps) and try again. The result is a list of steps. A separate program, {L('chain_verify.py')}, which shares no code with the search, checks the conditions of the three theorems on every step.

Random primes, {CH['big'][0]['primes']} of each size up to 60 digits and {[x for x in CH['big'] if x['digits'] == max(digs)][0]['primes']} above ({L('chain_test_big.py')} → {L('chain_test_big.txt')}, {L('chain_test_big2.log')}); the checker accepted every chain that was built, and PARI's isprime, which is a proof, agreed on every prime up to {max(int(k) for k in CH['pari'])} digits. Times are for plain Python.

| Digits | {' | '.join(n for _, n in modes)} | Median steps (first column) | Median time, s |
|---|---|---|---|---|---|---|
{rows}

On every number from 2 to {f(sm['every size of wheel']['top'])}, with trial division allowed only below 50, the chain proved all {f(sm['every size of wheel']['proved'])} primes and no other number; with curves only it proved {f(sm['curves only']['proved'])} and found no chain for {sm['curves only']['undecided']} prime ({sm['curves only']['und'][0]}) ({L('chain_test.py')} → {L('chain_test_small.txt')}). The page's code, run under node, proved all {CH['js_proved']} random primes of 20 to {CH['js_maxdigits']} digits given to it (median {j[100]['median_ms']} ms at 100 digits, {j[CH['js_maxdigits']]['median_ms']} ms at {CH['js_maxdigits']}) and all primes up to {f(jr[0]['numbers'] + 1)}, and the Python checker accepted every one of its chains ({L('chain_js_check.py')} → {L('chain_js_check.txt')}; the primes are in the files primes_N.txt).

Numbers that are not prime ({L('chain_adversary.py')} → {L('chain_adversary.txt')}): all {adv['spsp']['count']} strong pseudoprimes to base 2 below {f(adv['spsp']['top'])} enter the search, and for every one a turn failed that every prime makes or a factor appeared; of {npq} composites p(2p − 1) of 20 to 60 digits that pass the doubling wheel none got a chain; of {f(adv['mut']['made'])} tampered certificates (one number changed, or the top number replaced by a composite) the checker accepted {adv['mut']['accepted']}, all for a number that is prime (a changed base can still be a valid base).

What this is not. It is not new: it is elliptic curve primality proving, with few curve families and plain arithmetic. It is not one wheel but a chain, and finding the chain is a search for which no bound is proved; with these {CH['families']} families the table shows where it starts to fail. The page stops at 150 digits for the curve steps. The method itself has no such limit: the largest proof of this kind that we found reported has 109,297 digits.

'''
before('**What this is and is not.**', text)
pc = open('page_check_17.txt').read()
assert 'machine 17: Lamps lit so far: 669. Mistakes so far: 0. script errors 0' in pc
rep("and with 5,000 wheels the page's own code lights 669 lamps and makes no mistake on each of machines 10 to 16 (", "machine 17 is the chain of wheels (" + L('chain_summary.py') + " → " + L('chain.txt') + ", " + L('page_check_17.txt') + "); and with 5,000 wheels the page's own code lights 669 lamps and makes no mistake on each of machines 10 to 17 (")
open('section20.md', 'w', encoding='utf-8').write(s)
print(len(s))
