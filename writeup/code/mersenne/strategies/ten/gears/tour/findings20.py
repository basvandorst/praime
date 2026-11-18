"""Write the FINDINGS.txt block for section 20 from tournament.json."""
import json
T=json.load(open('tournament.json')); f=lambda n:format(n,',')
import re
FC=re.search(r'ring of 5 cells: guess fails for (\d+) of (\d+)',open('factor_classes.txt').read()).groups()
B64=re.findall(r': (\d+) numbers from ([0-9.e+]+) to ([0-9.e+]+); home on Pascal\'s wheel: (\d+); home on the two-cell wheel: (\d+)',open('beyond64.txt').read())
rs=sorted(T['census'],key=int); S=T['single']
blk=f'''
From the contest of wheels that decide alone (18 November 2026; writeup/code/mersenne/strategies/ten/gears/tour)
- One evaluation makes a search finite. A ring of r cells in which every cell is added to its neighbour n times, mod n, holds (x + 1)^n mod (n, x^r - 1) and reads x^n + 1 for every prime n. At x = 1 (the sum of the cells) this is 2^n = 2 (mod n), so a composite that passes is a base-2 pseudoprime, and the published list of those below 2^64 ({f(T['list_count'])} numbers) holds every candidate. Look for the evaluation that maps a new test onto a list that exists.
- A number cheats only on a wheel it fits. All {f(T['fit_liars'])} passes (number, ring) on rings of 2 to 40 cells below 1e10 have n^2 = 1 (mod r); on 5 cells {T['ring5_classes_1e10']['1']} have n = 1 and {T['ring5_classes_1e10']['4']} have n = -1 (mod 5). With r the smallest ring that n does not fit (gcd(n, r) = 1, n^2 != 1 mod r; 5 to {T['misfit']['max_ring']} cells) no candidate below 2^64 passes, and none passes any ring of {", ".join(rs)} cells that it does not fit. This is Agrawal's conjecture checked to 2^64 for every r up to {T['all_r_upto']}, given the list (the checks we found in the literature: 1e10 for r < 100, 1e11 for r = 5, a search to 1e17).
- Two cells are enough so far: (1 + sqrt c)^n = 1 - sqrt c (mod n), c the smallest odd prime with (c/n) = -1 (Khashin's Frobenius test), has no composite solution among the {f(S['duo']['tested1e9'])} odd composites below 1e9 and none among the base-2 pseudoprimes below 2^64. A fixed wheel (one polynomial for every n) always failed in our runs or is expected to: it fits some numbers. Above 2^64, {f(sum(int(b[0]) for b in B64))} numbers built to fool the doubling wheel (Chernick's Carmichael numbers and p(k(p - 1) + 1), up to 6e24) pass neither kept wheel.
- Prime squares: (x + 1)^(p^2) = x^(p^2) + 1 + p h(x^p) (mod p^2) with p h(x) = (x + 1)^p - x^p - 1, so p^2 passes on r cells exactly when h = 0 mod (p, x^r - 1). h(1) = (2^p - 2)/p forces a Wieferich prime; x(x + 1)(x^2 + x + 1) divides h; for 1093 and 3511 only r = 2, 3, 4 pass (r <= 200).
- A prediction that failed: "every prime factor p of a passing n has a power p^j = n (mod r)" is false for {FC[0]} of the {FC[1]} passes on 5 cells (a factor p = 1 with n = -1 mod 5: zeta^-1 + 1 = (zeta + 1)/zeta makes n = -1 as easy as n = 1).
'''
open('findings20.txt','w').write(blk); print(blk)
