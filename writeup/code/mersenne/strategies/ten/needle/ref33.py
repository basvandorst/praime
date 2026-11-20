"""Independent reference for prime33.html: Fermat quotients q(p) = (2^(p-1) - 1)/p mod p for odd primes p <= 10^6 with
plain Python pow (no gmpy2), their sum, the Wieferich primes, the count with min(q, p - q)/p < 1e-3 and 20 sector counts;
compared with the page's own numbers (page33.json).  usage: ref33.py page33.json out.txt"""
import sys, json
from sympy import primerange
B = json.load(open(sys.argv[1])); out = open(sys.argv[2], 'w')
n = 0; sm = 0; wf = []; near = 0; bins = [0] * 20
for p in primerange(3, 10 ** 6 + 1):
    q = (pow(2, p - 1, p * p) - 1) // p; n += 1; sm += q
    if q == 0: wf.append(p)
    bins[min(19, q * 20 // p)] += 1
    if min(q, p - q) * 1000 < p: near += 1
ok = (n, sm, wf, near, bins) == (B['n'], B['sum'], B['wief'], B['near'], B['bins'])
s = 'odd primes to 10^6: %d, sum of q(p) %d, Wieferich %s, near %d, sectors %s; page agrees: %s' % (n, sm, wf, near, bins, ok)
print(s); out.write(s + '\n')
