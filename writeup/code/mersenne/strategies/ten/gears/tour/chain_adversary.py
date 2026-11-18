import json, random, sys, time
import gmpy2
import cypari2
from chain_wheel import Prover, Composite, sprp
from chain_verify import verify

pari = cypari2.Pari()
pari.allocatemem(200 * 10**6)
rng = random.Random(77)
lines = []
out = {}


def say(s):
    print(s, flush=True)
    lines.append(s)


def attempt(n, **kw):
    P = Prover(seed=n % 1000003, **kw)
    try:
        c = P.prove(n)
    except Composite:
        return 'composite'
    return c


TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 50000000
t = time.time()
sp = [n for n in range(9, TOP, 2) if gmpy2.is_strong_prp(n, 2) and not gmpy2.is_prime(n)]
res = {'composite': 0, 'none': 0, 'chain': 0}
for n in sp:
    c = attempt(n, stop=50)
    res['chain' if isinstance(c, list) else 'composite' if c == 'composite' else 'none'] += 1
say('strong pseudoprimes to base 2 below %d (they pass the first screen): %d; the search found a factor or a failed turn for %d, found no chain for %d, built a chain for %d (%.0f s)' % (TOP, len(sp), res['composite'], res['none'], res['chain'], time.time() - t))
out['spsp'] = dict(top=TOP, count=len(sp), **res)

for d in (20, 30, 40, 60):
    found = []
    tries = 0
    while len(found) < 20:
        p = int(gmpy2.next_prime(rng.randrange(10**(d // 2 - 1), 10**(d // 2))))
        tries += 1
        q = 2 * p - 1
        if gmpy2.is_prime(q) and sprp(p * q):
            found.append(p * q)
    res = {'composite': 0, 'none': 0, 'chain': 0}
    for n in found:
        c = attempt(n)
        res['chain' if isinstance(c, list) else 'composite' if c == 'composite' else 'none'] += 1
    say('composites p(2p - 1) of about %d digits that pass the doubling wheel: %d; factor or failed turn %d, no chain %d, chain built %d' % (d, len(found), res['composite'], res['none'], res['chain']))
    out['pq%d' % d] = res

sq = [1093**2, 3511**2]
say('squares of the two Wieferich primes: ' + ', '.join('%d -> %s' % (n, 'chain' if isinstance(attempt(n, stop=50), list) else 'no chain') for n in sq))

mut = acc = acc_bad = same = 0
for d in (25, 40, 60):
    for _ in range(12):
        p = int(gmpy2.next_prime(rng.randrange(10**(d - 1), 10**d)))
        cert = attempt(p)
        if not isinstance(cert, list):
            continue
        assert verify(p, cert)
        for i, st in enumerate(cert):
            for k in st:
                if k in ('t', 'D'):
                    continue
                for delta in (1, -1, 2, rng.randrange(3, 10**6)):
                    c2 = json.loads(json.dumps(cert))
                    if k == 'f':
                        j = rng.randrange(len(st['f']))
                        c2[i]['f'][j][rng.randrange(2)] += delta
                    else:
                        c2[i][k] += delta
                    top = c2[0]['N']
                    mut += 1
                    if verify(top, c2):
                        acc += 1
                        if int(pari.isprime(top)) != 1:
                            acc_bad += 1
        c3 = json.loads(json.dumps(cert))
        for comp in (p + 2, p * 3, p * p):
            if int(pari.isprime(comp)) == 1:
                continue
            c3[0]['N'] = comp
            mut += 1
            if verify(comp, c3):
                acc += 1
                acc_bad += 1
say('tampered certificates: %d made (one number changed in one step, or the top number replaced by a composite); the checker still accepted %d, and of those %d were for a number that is not prime' % (mut, acc, acc_bad))
out['mut'] = dict(made=mut, accepted=acc, accepted_not_prime=acc_bad)
json.dump(out, open('chain_adversary.json', 'w'))
open('chain_adversary.txt', 'w').write('\n'.join(lines) + '\n')
