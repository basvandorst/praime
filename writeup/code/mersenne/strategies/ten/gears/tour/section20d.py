"""Fourth round of section 20: adds the test with a receipt (machine 18) to section20.md.
Every number comes from receipt.json, which receipt_summary.py builds from the output files of the test scripts."""
import json
RC = json.load(open('receipt.json'))
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


rows = '\n'.join('| %s | %s | %d | %s | %.2f | %.2f | %s and %d | %.3f | %.0f |' % (f(r['p']), r['group'], r['k'], f(r['squarings']), r['test_s'], r['receipt_s'], f(r['check_squarings']), r.get('small_powers', 2 * r['k']), r['check_s'], r['ratio']) for r in RC['timing'])
a = RC['agree']
t = RC['tamper']
fg = RC['forgery']
sf = RC['safe']
f3 = [d for d in fg if 'order 3' in d['factor']][0]
f4 = [d for d in fg if 'order 4' in d['factor']][0]
frows = '\n'.join('| %s | %s | %d | %d | %d | %d |' % ('one cell' if d['group'] == 'one' else 'two cells', d['factor'], d['exps'], d['tried'], d['accepted'], d['changed']) for d in fg)
text = f'''**A test with a receipt (the halving protocol and its weakness are literature; the runs are computed).** For a Mersenne number the Lucas–Lehmer test cannot be beaten by any wheel of this section; what can be saved is the second run that a result needs before it is trusted. The statement to check is y = x^(2^T) in a group: T squarings.

The protocol (Pietrzak, 2018). If T is odd, replace y by y² and T by T + 1. The tester sends the middle value μ = x^(2^(T/2)). Then y = μ^(2^(T/2)), and for every number r: (x^r · μ)^(2^(T/2)) = μ^r · y, which is again a statement of the same kind with half the turns. The number r is taken from a hash of x, y, μ and T (64 bits of SHA-256 here), so the tester cannot choose it. After k halvings the checker makes the remaining T/2^k squarings itself; in all it makes T/2^k squarings and 2k powers with exponents of 64 bits. The tester needs no second pass. With h_1, …, h_k the lengths of the halves, the i-th middle value is f(i − 1, h_i), where f(0, m) = x^(2^m) and f(j, m) = f(j − 1, m)^(r_j) · f(j − 1, m + h_j); so it keeps the residues at the 2^k − 1 positions h_i + (a sum of some of h_1, …, h_(i−1)) and builds the receipt from them.

Two groups ({L('receipt.py')}). One cell: x = 3 mod n = 2^p − 1 and T = p; for a prime n, 3^(n+1) = 3^(2^p) ≡ 9, so y = 9 is a probable-prime test (not a proof). Two cells: x = 2 + √3, which has norm 1, and T = p − 2; squaring is (a, b) → (2a² − 1, 2ab), because a² − 3b² = 1. Twice the first cell is x^(2^i) + x^(−2^i), the Lucas–Lehmer sequence 4, 14, 194, …, so the first cell is 0 after p − 2 squarings exactly when 2^p − 1 is prime (the Lucas–Lehmer theorem, p odd). This is machine 11's wheel with step 3 again, and its verdict is a proof.

Agreement ({L('receipt_test.py')} → {L('receipt_test.txt')}): for all {a['primes']} odd primes p up to {f(a['top'])} both forms give the verdict of the plain Lucas–Lehmer test, and every honest receipt is accepted; the exponents found are the 19 odd Mersenne exponents below 5,000. For no composite exponent below 2,000 does 2^n − 1 pass the one-cell test.

Cost ({L('receipt_timing.py')} → {L('receipt_timing.txt')}; Python with gmpy2, the smaller time of two runs on a shared machine, and for the test the smallest of its four runs; for each p the receipt of 6 or 8 numbers that gave the faster check). The counts of squarings are exact; the times are one machine's.

| p | Form | Halvings | Turns of the test | Test, s | Receipt, s | Check: squarings and small powers | Check, s | Test / check |
|---|---|---|---|---|---|---|---|---|
{rows}

For small p the 2k small powers cost more than the turns that are saved, so the gain grows with p. The project that searches for Mersenne primes reports that its receipts are checked with less than 0.5% of the work of a test (literature, mersenne.org); it uses the one-cell form.

Tampering, p = 9,689, receipts of 6 numbers, each form: a wrong result was accepted in {t['two']['wrong_y_accepted']} of {t['two']['wrong_y_tried']} tries (the true result plus 1, and 50 random values), a changed receipt number in {t['two']['wrong_mu_accepted']} of {t['two']['wrong_mu_tried']}, and after a single bit was flipped at a random turn of the test the receipt was accepted in {t['two']['flips_accepted']} of {t['two']['flips']} runs (the result was wrong in all of them).

A forged receipt (the weakness is literature: Hoffmann, Hubáček, Kamath and Pietrzak, 2023, who state it for Proth numbers and repair it; the runs are ours). Let ε be an element with ε^d = 1 for a small d, and suppose the tester knows the true y and μ. It reports y·ε and sends μ·δ. The check of the first halving, (x^r · μδ)^(2^h) = (μδ)^r · yε with h = T/2, holds exactly when δ^(2^h − r) = ε. Taking δ = ε^t, this is t·(2^h − r) ≡ 1 (mod d), where r depends on t through the hash; each t without a common factor with d can be tried, and each fits with a chance of about 1/d, because r is a hash value (heuristic). After that the new statement is true and the rest of the receipt is honest. If T is odd the first step squares y, so ε is replaced by ε² and the factor −1 disappears at once: the receipt cannot see the sign. For a cube root of 1 the heuristic gives 1 − (2/3)² = 56% and we saw {100 * f3['accepted'] / f3['tried']:.0f}%; for the factor of order 4 on two cells (ε² = −1, one try, chance 1/2) we saw {100 * f4['accepted'] / f4['tried']:.0f}%. In the table an attempt is one exponent with one label added to the hash (40 labels per exponent, to count how often the hash falls right; a real receipt has no such label, so a forger has one try per factor).

| Form | Factor ε | Prime 2^p − 1 used | Attempts | Accepted by the checker | Verdict changed from prime to not prime |
|---|---|---|---|---|---|
{frows}

On two cells the factor −1 leaves the first cell at 0, so the verdict stays; the factor of order 4 does not. So a receipt of this kind shows that the turns were made by someone who knows the true result. It does not show that the reported result is the true one.

A guard (ours; tested, not proved). Compare not the result but a power of it that removes every factor of small order: on one cell y^L against 9^L, with L = 2 times the prime factors below 2^16 of (n − 1)/2; on two cells whether y^(2^{sf['S']}) = 1. A result that passes this weaker comparison is sent to a full test. All {sf['one_forged'] + sf['two_forged']} forged results tried were caught this way ({sf['one_forged_flagged']} of {sf['one_forged']} on one cell, {sf['two_forged_flagged']} of {sf['two_forged']} on two), and of the {sf['composites']} composite 2^p − 1 with prime p < 2,000 none was sent to a full test without need.

The page's code and the Python code were checked against each other ({L('receipt_js_check.py')} → {L('receipt_js_check.txt')}): for {RC['js']['exps']} exponents from 3 to 9,689 the page's receipt is identical to Python's and accepted by it, and the page's checker accepts {RC['js']['honest']} honest Python receipts and rejects {RC['js']['tampered']} tampered ones. The page's own SHA-256 agrees with node's on 6 strings.

'''
before('**What this is and is not.**', text)
assert 'machine 18: Lamps lit so far: 20. Mistakes so far: 0. script errors 0' in open('page_check_18.txt').read()
rep("and with 5,000 wheels the page's own code lights 669 lamps and makes no mistake on each of machines 10 to 17 (", "machine 18 is the two-cell wheel with a receipt, where wheel p decides 2^p − 1 (" + L('receipt_summary.py') + " → " + L('receipt.txt') + "; with 5,000 wheels it lights the 20 Mersenne exponents below 5,000 and nothing else, " + L('page_check_18.txt') + "); and with 5,000 wheels the page's own code lights 669 lamps and makes no mistake on each of machines 10 to 17 (")
open('section20.md', 'w', encoding='utf-8').write(s)
print(len(s))
