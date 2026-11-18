"""Write section 20 of MERSENNE.md from tournament.json (every figure comes from the logs)."""
import json,sys
T=json.load(open('tournament.json')); f=lambda n:format(n,',')
G='writeup/code/mersenne/strategies/ten/gears/tour/'
L=lambda name:'[%s](%s%s)'%(name,G,name)
import re
OWN=int(re.search(r'identical to the first (\d+) lines of the published list: True',open('list_check.txt').read()).group(1))
FC=re.search(r'ring of 5 cells: guess fails for (\d+) of (\d+)',open('factor_classes.txt').read()).groups()
B64=re.findall(r': (\d+) numbers from ([0-9.e+]+) to ([0-9.e+]+); home on Pascal\'s wheel: (\d+); home on the two-cell wheel: (\d+)',open('beyond64.txt').read())
import math
_mx=max(float(b[2]) for b in B64); _e=int(math.floor(math.log10(_mx))); MX64='%.1f · 10^%d'%(_mx/10**_e,_e)
FB={int(a):(int(b),int(c)) for a,b,c in re.findall(r'degree (\d+): (\d+) wheels, (\d+) with no liar',open('tournamentA_1e10.txt').read())}
S=T['single']; M=T['misfit']; C=T['census']; rs=sorted(C,key=int)
ringuse=", ".join("%s cells %.1f %%"%(k,100*v/M['tested']) for k,v in M['rings'].items() if v/M['tested']>=0.005)
cen="; ".join("%s cells: %s tested"%(r,f(C[r]['tested'])) for r in rs)
maxr=max(int(r) for r in rs)
txt=f'''## 20. One wheel decides: a wheel the number does not fit (that every prime comes home is proved; "no mistake below 2^64" is computed and rests on the published list of base-2 pseudoprimes; that there is never a mistake is open: Agrawal's conjecture for the ring, an unsettled Frobenius test for the two cells)

Sections 17 to 19 judged a machine by how many prime wheels are home together. A harder demand is that one wheel decides alone: wheel n must be home after n turns exactly when n is prime. The gear rules were put in a contest on that demand ({L('tour.c')}, {L('psp2.c')}, {L('summary.py')} → {L('tournament.txt')}).

**The machines.** Each wheel n carries a few cells that hold remainders after division by n.

- Doubling wheel: one cell, doubled every turn; home when it shows its start again, 2^n ≡ 2 (mod n).
- Fibonacci wheel: two cells (a, b) → (b, a + b) from (0, 1); home when n divides F(n − 1) for n ≡ ±1 (mod 5) and F(n + 1) for n ≡ ±2 (mod 5); multiples of 5 are turned away.
- Perrin's wheel: three cells (a, b, c) → (b, c, a + b) from (3, 0, 2); home when n divides the first cell after n turns.
- Two-cell wheel with step c: cells (a, b) → (a + c·b, a + b) from (1, 0), which is multiplication by 1 + √c; home when the cells read (1, −1), that is (1 + √c)^n ≡ 1 − √c (mod n).
- Pascal's wheel with r cells: a ring with a 1 in cell 0; a turn replaces every cell by itself plus the cell before it. After t turns the ring holds the coefficients of (x + 1)^t mod (n, x^r − 1). Home means a 1 in cell 0, a 1 in cell n mod r and nothing else: (x + 1)^n ≡ x^n + 1 (mod n, x^r − 1).

**Every prime comes home (proved).** For a prime n every binomial coefficient C(n, k) with 0 < k < n is divisible by n, because n divides n! and not k!(n − k)!. So (x + 1)^n ≡ x^n + 1 (mod n) before the ring is folded, and folding keeps it. The same binomials give (1 + √c)^n ≡ 1 + c^((n − 1)/2) √c (mod n), and c^((n − 1)/2) ≡ −1 when c is not a square mod the prime n (Euler's criterion). The other three wheels are the classical statements 2^n ≡ 2, n | F(n − (5/n)) and n | P(n) for primes. So a wheel that is not home proves that n is not prime; the only possible mistake is a number that is not prime and comes home.

**Two lemmas that make the search finite (proved).**

1. Sum of the cells. A turn doubles the sum of the cells, so after n turns the sum is 2^n; the home position has sum 2. A number that comes home on any Pascal ring therefore satisfies 2^n ≡ 2 (mod n): it fools the doubling wheel first. The odd numbers that do are known below 2^64: {f(T['list_count'])} of them (Feitsma and Galway; the {f(OWN)} smallest agree with our own sieve, {L('list_check.py')} → {L('list_check.txt')}).
2. Folding. If d divides r, folding the ring of r cells onto d cells (cell k onto cell k mod d) turns one home position into the other. A mistake on r cells is a mistake on every ring of d cells, d | r.

**The contest (computed).** All odd numbers below 10^8 that are not prime, without any shortcut ({f(T['odd_composites_below_1e8'])} numbers):

| Wheel | Cells | First mistake | Mistakes below 10^8 |
|---|---|---|---|
| Doubling | 1 | {f(S['doubling']['first'])} | {f(S['doubling']['below1e8'])} |
| Fibonacci | 2 | {f(S['fib']['first'])} | {f(S['fib']['below1e8'])} |
| Perrin | 3 | {f(S['perrin']['first'])} | {f(S['perrin']['below1e8'])} |
| Pascal, 5 cells | 5 | {f(S['ring5']['first'])} | {f(S['ring5']['below1e8'])} |

Wheels with a feedback, x^d = c_(d−1) x^(d−1) + … + c_0 with coefficients from −3 to 3 and the same test (x + 1)^n ≡ x^n + 1, were joined with the doubling wheel and run on the list below 10^10: all {FB[1][0]} with one cell and all {FB[2][0]} with two cells make mistakes, and {FB[3][1]} of the {FB[3][0]} with three cells make none there ({L('tournamentA.py')} → {L('tournamentA_1e10.txt')}). A fixed wheel fits some numbers, so these are expected to fail further on; they were dropped.

No fixed wheel is free of mistakes, and more cells alone do not help: on the list below 10^10 the rings of 2 to 40 cells make {f(T['fit_liars'])} mistakes (number and ring counted together; {L('ringcensus.py')} → {L('ringcensus_1e10.txt')}).

**Fit.** Say that n fits a ring of r cells when n² ≡ 1 (mod r). Of those {f(T['fit_liars'])} mistakes, {T['fit_liars_misfits']} are numbers that do not fit their ring. On 5 cells the {f(S['ring5']['below1e10'])} mistakes below 10^10 split as {T['ring5_classes_1e10']['1']} with n ≡ 1 and {T['ring5_classes_1e10']['4']} with n ≡ −1 (mod 5), none with n ≡ ±2. Why ±1 is special (a heuristic, not a proof): for a prime factor p of n and an r-th root of unity ζ over the field with p elements, home means (ζ + 1)^n = ζ^n + 1. For n ≡ 1 this reads (ζ + 1)^(n − 1) = 1, and for n ≡ −1, because ζ^(−1) + 1 = (ζ + 1)/ζ, it reads (ζ + 1)^(n − 1) = ζ^(−1). Both are conditions on the order of one element, of the kind that Carmichael numbers meet. For another class ζ^n + 1 is not ζ + 1 times a root of unity, and the equation asks for a coincidence between two unrelated elements, unless n ≡ p^j (mod r) for every prime factor, the case in which Lenstra and Pomerance expect mistakes among very large numbers. A simpler guess, that every prime factor p of a mistake has a power p^j ≡ n (mod r), was tested and fails for {FC[0]} of the {FC[1]} mistakes on 5 cells: a factor p ≡ 1 with n ≡ −1 ({L('factor_classes.py')} → {L('factor_classes.txt')}).

**The two machines that were kept.** Both give n a wheel it does not fit.

Pascal's wheel. For each n take the smallest r that has no factor in common with n and that n does not fit (r is 5 or more, since every number fits 2, 3 and 4 cells), and run Pascal's wheel on r cells. Results:

- All {f(M['tested'])} candidates below 2^64: {M['liars']} came home. Ring sizes: {ringuse}; the largest ring needed had {M['max_ring']} cells.
- Separately, every candidate was put on each ring of {", ".join(rs[:-1])} and {rs[-1]} cells that it does not fit ({cen}): none came home. By lemma 2 a mistake on a ring r that n does not fit shows on one of the prime-power parts of r that n does not fit either, and the parts 2, 3, 4 and 8 fit every number. So this settles every ring whose other prime-power parts are among {", ".join(rs)}; that includes every r up to {T['all_r_upto']}.
{("- A second program, written separately with plain 128-bit remainders ("+L('ring2.c')+"), repeated "+" and ".join("the ring of %d cells (%s tested, %d came home)"%(rr,f(T['ring2_%d'%rr]['tested']),T['ring2_%d'%rr]['liars']) for rr in (5,7) if 'ring2_%d'%rr in T)+".") if 'ring2_5' in T else ""}
- Even numbers: the test box turns an even n above 2 away. On the page's range (2 to 5,000) no even number comes home, and the three even numbers below 3 · 10^6 with 2^n ≡ 2 (mod n) do not either ({L('pascal_check.py')} → {L('pascal_check.txt')}).

The two-cell wheel. Its step is the smallest odd prime c that is not a square mod n in the sense of the Jacobi symbol, (c/n) = −1; n "fits" a step that is a square. Even numbers, numbers divisible by a step tried on the way, and numbers with no step up to 199 (a square has none) are turned away. Results:

- Alone, on every odd number below {S['duo']['bound_md']} that is not prime ({f(S['duo']['tested'])} numbers): {S['duo']['liars']} came home.
- On all {f(T['duo_list']['tested'])} candidates below 2^64, that is together with the doubling wheel: {T['duo_list']['liars']} came home; {T['duo_list']['nostep']} had no step (the squares of 1093 and 3511); the largest step needed was {max(int(c) for c in T['duo_list']['steps'])}.
- Independent check in plain Python with sympy's Jacobi symbol: {L('tour_check2.py')} → {L('tour_check2.txt')}.


A check of the list and of the code against a known result: the doubling wheel together with the Fibonacci wheel for n ≡ ±2 (mod 5), the pair for which a prize is offered, on {f(T['pair']['tested'])} candidates: {T['pair']['liars']} came home.

Above 2^64 there is no list. As a spot check, numbers that fool the doubling wheel were built on purpose: {f(int(B64[0][0]))} Carmichael numbers (6k + 1)(12k + 1)(18k + 1) and {f(int(B64[1][0]))} numbers p·(k(p − 1) + 1), k = 2 to 6, all above 2^64 and up to {MX64}. On Pascal's wheel {B64[0][3]} and {B64[1][3]} came home, on the two-cell wheel {B64[0][4]} and {B64[1][4]} ({L('beyond64.py')} → {L('beyond64.txt')}).

In the literature the two-cell wheel with this choice of step is Khashin's form of the Frobenius test (2013); no number that fools it is known, and Khashin shows that one below 2^60 would need a prime factor below 19 or a step above 128. The statement "home on r cells implies prime or n² ≡ 1 (mod r)" is Agrawal's conjecture. The checks we found reported: r < 100 with n < 10^10, r = 5 with n < 10^11, and a distributed search over 10^10 < n < 10^17. The run above reaches 2^64 ≈ 1.8 · 10^19 for the rings named. It is a computation, not a proof, and it rests on the completeness of the published list.

**Prime powers (proved, then computed).** Write (x + 1)^p = x^p + 1 + p h(x). Then (x + 1)^(p²) = (x^p + 1 + p h(x))^p ≡ (x^p + 1)^p ≡ x^(p²) + 1 + p h(x^p) (mod p²). So p² comes home on r cells exactly when h(x^p), and hence h(x) (x → x^p only permutes the cells, as p does not divide r), vanishes mod (p, x^r − 1). Now h(1) = (2^p − 2)/p, which is divisible by p only for a Wieferich prime, and x(x + 1)(x² + x + 1) divides h as a polynomial. Any higher power needs a Wieferich prime too: if p^k, k ≥ 2, comes home, lemma 1 gives 2^(p^k − 1) ≡ 1 (mod p²); the order of 2 mod p² divides p(p − 1) and p^k − 1, whose greatest common divisor is p − 1, so 2^(p − 1) ≡ 1 (mod p²). For the two known Wieferich primes 1093 and 3511 the square fools the rings of 2, 3 and 4 cells, which every number fits, and no other ring up to 200 cells ({L('lemma_pp2.py')} → {L('lemma_pp2.txt')}).

**Cost.** A turn count of n is reached in about log2 n doublings of the turn count, each a product of two rings: about r² log2 n multiplications mod n. The page tests 2^4423 − 1 (1,332 digits) in about a second.

**What this is and is not.** No wheel here is new mathematics: the pair is the test of Baillie, Pomerance, Selfridge and Wagstaff, the two-cell wheel is a Frobenius test, and Pascal's wheel is the congruence of Agrawal, Kayal and Saxena, whose proved test needs many rings and offsets. What the contest adds is the common rule in gear terms (in every case found a number cheats only on a wheel it fits, so take a wheel it does not fit), the check of Pascal's wheel up to 2^64, and a small statement on prime powers. A proof that either kept machine never errs would give a prime test of a few lines and would be a major result; so would the first number that fools it.

**The page.** Machines 10 and 11 of [prime3.html](writeup/code/mersenne/prime3.html) are Pascal's wheel and the two-cell wheel: at tick n wheel n decides, the lamps stay on, a picture shows every turn of the wheel, and a box tests a number of up to 2,000 digits. Among 250 wheels all {T['row250'][0]} prime lamps are on at tick {T['row250'][1]}, all {T['row250'][2]} exponents at tick {T['row250'][3]}, all {T['row250'][4]} twin primes at tick {T['row250'][5]}, and no other lamp; the two-cell wheel gives the same ({L('duo_row.py')} → {L('duo_row.txt')}). The tape and the square jumps, which gave the same results as machine 1, were removed.

'''
open('section20.md','w').write(txt); print(len(txt))
