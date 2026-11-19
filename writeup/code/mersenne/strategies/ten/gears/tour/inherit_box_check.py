"""The box of machine 23 (page code, taken out of prime3.html by inherit_node.js) against the rule from the top down with PARI.
usage: inherit_box_check.py make numbers.txt   |   inherit_box_check.py check inherit_js.json"""
import sys, json, random
from inherit_big import rounds, pari

if sys.argv[1] == 'make':
    rng = random.Random(23)
    L = [int(pari.nextprime(k)) for k in range(3, 3000) if int(pari.nextprime(k)) == k]
    for digits, count in ((9, 150), (12, 150), (15, 100), (18, 100), (22, 40), (26, 30), (30, 20)):
        L += [int(pari.nextprime(rng.randrange(10 ** (digits - 1), 10 ** digits))) for _ in range(count)]
    L += [21767, 1688173, 2223990533, 3053740441, 3859111517, 944569134151125699012274669297]
    L += [76175908868713, 76790843706589, 299410817655227, 318142460798927, 38087954434357, 3173996202863]
    open(sys.argv[2], 'w').write('\n'.join(str(x) for x in L) + '\n')
else:
    J = json.load(open(sys.argv[2]))
    B = J['box']
    pd = []
    bad = und = part = 0
    for t, r, ok, p in B:
        if r < 0:
            und += 1
        elif p:
            part += 1
            pd.append(len(t))
            if r < rounds(int(t)) or not ok:
                bad += 1
        elif r != rounds(int(t)) or not ok:
            bad += 1
    print('box of machine 23: %d primes (all below 3,000; random ones of 9 to 30 digits; the first primes of 4, 5 and 6 rounds; four primes of 7 rounds): rounds that differ from the rule worked out with PARI, or a proof that failed: %d; primes where the page could not take a neighbour apart and says so: %d; not finished by the page: %d' % (len(B), bad, part, und))
    print('the primes where the page says so: they have %d to %d digits' % (min(pd), max(pd)))
    print('helpers of 2223990533 on the page, largest first: ' + ', '.join(str(x) for x in J['chain']))
