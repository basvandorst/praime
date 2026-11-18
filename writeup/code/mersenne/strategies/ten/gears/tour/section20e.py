"""Fifth round of section 20: adds the two-cell wheel from a start of norm 2 (machine 19) to section20.md.
Every number comes from norm2.json, which norm2_summary.py builds from the output files of the programs."""
import json
N2 = json.load(open('norm2.json'))
f = lambda n: format(n, ',')
G = 'writeup/code/mersenne/strategies/ten/gears/tour/'
L = lambda n: '[%s](%s%s)' % (n, G, n)
s = open('section20.md', encoding='utf-8').read()


def before(anchor, text):
    global s
    assert s.count(anchor) == 1, anchor[:60]
    s = s.replace(anchor, text + anchor)


def rep(old, new):
    global s
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new)


import re
_b = re.findall(r': (\d+) numbers from .*?home on the norm-2 wheel: (\d+)', open('norm2_beyond.txt').read())
BY = (sum(int(x) for x, y in _b), sum(int(y) for x, y in _b))
assert len(_b) == 3
JS = re.search(r'(\d+) numbers \(.*?different: (\d+)', open('norm2_js_check.txt').read()).groups()
assert JS[1] == '0'
li = N2['list']
SEC = (" A second program, written separately ({L('norm2_second.c')}: plain remainders instead of Montgomery form, the power taken from the lowest bit up), read the whole list again: on a wheel {f(li['second']['tested'])}, came home {li['second']['home']}, and every count equals the first run.".replace("{L('norm2_second.c')}", L('norm2_second.c')).replace("{f(li['second']['tested'])}", f(li['second']['tested'])).replace("{li['second']['home']}", str(li['second']['home']))) if 'second' in li else ''
rg = N2['range']
lm = N2['lemma']
bp = N2['bench']['primes']
bo = N2['bench']['odd']
k1 = 'doubling wheel'
k7 = 'seven strong tests'
kn = 'norm-2 wheel'
ko = "machine 11's wheel"
text = f'''**A two-cell wheel that is settled alone below 2^64 (the norm lemma is literature; the run is computed and rests on the published list).** Pascal's wheel is settled alone below 2^64 by lemma 1 (the sum of its cells is the doubling wheel), but it needs 5 to 47 cells. Machine 11 has two cells and was checked below 2^64 only together with the doubling wheel. A different start gives a wheel of two cells that needs no second test.

The wheel. Let n be odd and not a perfect square. For a = 1, 2, 3, … take c = a² − 2 (from a = 2 on) and then c = a² + 2. If c and n share a factor and n does not divide c, n is not prime. Stop at the first c with Jacobi symbol (c/n) = −1. The two cells start at (a, 1), that is a + √c, and n is called prime when (a + √c)^n ≡ a − √c (mod n). The starts are 1 + √3, 2 + √2, 2 + √6, 3 + √7, 3 + √11, …

1. Every prime comes home (proved). For a prime n and (c/n) = −1 the numbers x + y√c mod n form a field with n² elements, in which the n-th power is the map √c → −√c. A suitable a exists: for a prime n > 3 the sum of the Legendre symbols of a² + 2 over all a mod n is −1, so a² + 2 is a non-residue for at least (n − 1)/2 values of a; n = 3 uses 2 + √2.
2. The norm lemma (literature: a Frobenius pseudoprime for x² − Px + Q is a Fermat pseudoprime to base Q; here P = 2a, Q = ±2). The map √c → −√c is a ring automorphism mod n. Apply it to (a + √c)^n ≡ a − √c: also (a − √c)^n ≡ a + √c. Multiply the two: (a² − c)^n ≡ a² − c (mod n). Now a² − c = 2 or −2 and n is odd, so 2^n ≡ 2 and 2^(n−1) ≡ 1 (mod n). A composite n that comes home is a base-2 pseudoprime.
3. The run ({L('norm2.c')}; the list was streamed from the download and not stored, sha256 {li['sha256'][:8]}…{li['sha256'][-4:]} as before). All {f(li['lines'])} numbers of the list were read (last {f(li['last'])}). {f(li['tested'])} were put on their wheel and {li['home']} came home; {f(li['factor'])} share a factor with a small c and {len(li['squares'])} are perfect squares ({', '.join(f(x) for x in li['squares'])}), so they are turned away before any turn. The largest a was {li['amax']}. {L('norm2_list_check.py')} → {L('norm2_list_check.txt')} recomputes every 4,000th line in plain Python ({L('norm2.py')}): {f(li['samples'])} lines, {li['sample_bad']} differ.{SEC}

4. Why it behaves like a pair of wheels (proved). Write a² − c = 2s with s = 1 or −1. Then (a + √c)² = (a² + c) + 2a√c = 2ε with ε = (a² − s) + a√c, and ε has norm (a² − s)² − a²c = 1: the square of the start is twice a unit. So (a + √c)^n = (a + √c) · 2^((n−1)/2) · ε^((n−1)/2), and (a − √c)/(a + √c) = (a − √c)²/(2s) = s/ε. Home is therefore the single equation 2^((n−1)/2) · ε^((n+1)/2) ≡ s (mod n): Euler's test for 2 and a Lucas test for the unit ε, tied together. For the start 1 + √3 the unit is 2 + √3, the wheel of the Lucas–Lehmer test; for 2 + √2 it is 3 + 2√2, the square of 1 + √2. The pair "doubling wheel and a Lucas test" is the idea of Baillie, Pomerance, Selfridge and Wagstaff; here both turn on one shaft, and the first is not run separately.

So, given that the published list is complete: for every odd n below 2^64 this wheel alone says prime exactly when n is prime. The same lemma covers machine 11 when its step is 3, because 1 + √3 has norm −2; for the steps 5, 7, 11, … the norm is −4, −6, −10, … and the list does not apply.

Checks that do not use the list. The lemma was tested where it is not empty ({L('norm2_check.py')} → {L('norm2_check.txt')}): on wheels that n fits ((c/n) = +1), {lm['cases']} odd composites below {f(lm['top'])} satisfy (a + √c)^n ≡ a ± √c, and every one is a base-2 pseudoprime. And the wheel was put on every odd number from {rg['lo']} to 10^10: {f(rg['home'])} came home, and a sieve of Eratosthenes that shares nothing with the wheel ({L('pcount.c')}) counts {f(rg['primes'])} odd primes there. The largest a needed by a prime was {rg['amax']} (n = {f(rg['amaxn'])}).

The page's code for machine 19, run under node, gives the same start, the same two cells after n turns and the same verdict as {L('norm2.py')} on {f(int(JS[0]))} numbers (the sampled list lines, 3,000 primes and 3,000 random odd numbers of 19 digits, and every odd n up to 20,001): {JS[1]} differ ({L('norm2_js_check.py')} → {L('norm2_js_check.txt')}).

What is known. Khashin (2018) states that his own version of the two-cell test (steps −1 and 2 first, then 1 + √c) has no mistake below 2^64; his proof is a long computation on the possible prime factors and does not use the list. We did not find the route above (a start of norm 2, so that the list of base-2 pseudoprimes settles the wheel in one pass) written down; it is a direct consequence of the lemma, so it may well be known. Above 2^64 nothing is proved, and the reasons to expect a mistake of such tests at some very large number apply here too. As a small probe, {f(BY[0])} numbers above 2^64 that were built to fool the doubling wheel (Carmichael numbers of Chernick's form and pseudoprimes p(k(p − 1) + 1), up to 6 · 10^60) were put on the wheel: {BY[1]} came home ({L('norm2_beyond.py')} → {L('norm2_beyond.txt')}).

Cost ({L('norm2_bench.txt')}; one run in C, numbers between 2^62 and 2^64, time per number). On {f(bp[kn]['count'])} primes: one strong test to base 2 (the doubling wheel) {bp[k1]['ns']} ns, the norm-2 wheel {bp[kn]['ns']} ns, machine 11's wheel {bp[ko]['ns']} ns, the seven strong tests that are known to decide every number below 2^64 {bp[k7]['ns']} ns. On {f(bo[kn]['count'])} random odd numbers, where most are turned away early: {bo[k1]['ns']}, {bo[kn]['ns']}, {bo[ko]['ns']} and {bo[k7]['ns']} ns. The program is plain (four multiplications mod n for a squaring of the two cells); it is not tuned.

'''
before('**What this is and is not.**', text)
assert 'machine 19: Lamps lit so far: 669. Mistakes so far: 0. script errors 0' in open('page_check_19.txt').read()
rep("machine 18 is the two-cell wheel with a receipt, where wheel p decides 2^p − 1 (", "machine 19 is the two-cell wheel from a start of norm 2 (" + L('norm2_summary.py') + " → " + L('norm2.txt') + "; 669 lamps and no mistake with 5,000 wheels, " + L('page_check_19.txt') + "); machine 18 is the two-cell wheel with a receipt, where wheel p decides 2^p − 1 (")
rep("for the two cells it is an unsettled Frobenius test)", "for the two cells it is a Frobenius test, which a start of norm 2 settles below 2^64, given the list, and nothing settles above; the chain of wheels and the test with a receipt are literature and their runs are computed)")
open('section20.md', 'w', encoding='utf-8').write(s)
print(len(s))
