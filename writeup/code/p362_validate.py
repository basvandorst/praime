"""Second, independent check of puzzle 362 (plain integers, sympy).

Reads data/p362_terms.txt (n, a(n)) and checks, for every n:
  - a(n) is prime (sympy.isprime: a proof below 2^64, a PRP test above);
  - a(n) is the least prime q with  P(n-1) * (q - 1) > q * Q(n-1) / 2,
    where P/Q = prod_{i<n} (1 - 1/a(i)) kept as a pair of integers,
    i.e. the greedy rule prod (1 - 1/a(i)) > 1/2 at every step (A339344);
  - the page's nine terms are reproduced.
"""
from sympy import isprime, prevprime
page = [3, 5, 17, 257, 65537, 4294967311, 1229782942255939601,
        88962710886098567818446141338419231,
        255302062200114858892457591448999891874349780170241684791167583265041]
terms = [int(l.split()[1]) for l in open(__file__.rsplit('/code/', 1)[0] + '/data/p362_terms.txt')]
assert terms[:9] == page
P, Q = 1, 1                      # P/Q = prod (1 - 1/a(i)) so far, exact
for n, q in enumerate(terms, 1):
    assert isprime(q), n
    ok = lambda x: 2 * P * (x - 1) > x * Q          # P/Q * (1 - 1/x) > 1/2
    assert ok(q), n
    if n > 1:
        assert not ok(prevprime(q)), n              # the previous prime fails -> q is the least
    # the least prime > 2P/(2P - Q): every candidate between that bound and q is composite
    bound = (2 * P) // (2 * P - Q)                  # floor of the bound; q is the first prime above
    assert q > bound and all(not isprime(x) for x in range(bound + 1, q)), n
    P, Q = P * (q - 1), Q * q
    print(f"a({n}): {len(str(q))} digits ok")
print("the excess prod(1 - 1/a(i)) - 1/2 after a(12) is", f"{(2*P - Q)/(2*Q):.3e}" if False else "positive:", 2 * P > Q)
