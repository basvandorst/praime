import json, random, signal, sys, time
import gmpy2
import cypari2
from chain_wheel import Prover, Composite
from chain_verify import verify

pari = cypari2.Pari()
pari.allocatemem(200 * 10**6)
rng = random.Random(20261118)
sizes = [int(x) for x in sys.argv[1].split(',')]
count = int(sys.argv[2])
tlim = int(sys.argv[3])
modes = [('all sizes of wheel, 97 curve families', dict(hmax=4)),
         ('curves only, 97 families', dict(hmax=4, use_pm=False)),
         ('all sizes of wheel, 9 curve families', dict(hmax=1)),
         ('n - 1 and n + 1 only', dict(use_ec=False))]
if len(sys.argv) > 4:
    modes = [modes[int(i)] for i in sys.argv[4].split(',')]
tag = sys.argv[5] if len(sys.argv) > 5 else 'big'


class Timeout(Exception):
    pass


def alarm(*a):
    raise Timeout()


signal.signal(signal.SIGALRM, alarm)
rows = []
lines = []
for d in sizes:
    primes = []
    while len(primes) < count:
        p = int(gmpy2.next_prime(rng.randrange(10**(d - 1), 10**d)))
        if len(str(p)) == d:
            primes.append(p)
    pari_ok = sum(1 for p in primes if int(pari.isprime(p)) == 1) if d <= 120 else None
    for label, kw in modes:
        ok = fail = tout = bad = 0
        steps = []
        times = []
        nodes = []
        for p in primes:
            P = Prover(seed=p % 1000003, **kw)
            t = time.time()
            signal.alarm(tlim)
            try:
                cert = P.prove(p)
            except Timeout:
                cert = 'timeout'
            except Composite:
                cert = 'composite'
            signal.alarm(0)
            dt = time.time() - t
            if isinstance(cert, list):
                ok += 1
                if not verify(p, cert):
                    bad += 1
                steps.append(len(cert))
                times.append(dt)
                nodes.append(P.nodes)
            elif cert == 'timeout':
                tout += 1
            else:
                fail += 1
        med = lambda v: sorted(v)[len(v) // 2] if v else None
        row = dict(digits=d, mode=label, primes=count, proved=ok, no_chain=fail, timeout=tout, checker_rejects=bad, pari_agrees=pari_ok,
                   median_steps=med(steps), max_steps=max(steps) if steps else None, median_s=med(times), max_s=max(times) if times else None, median_nodes=med(nodes))
        rows.append(row)
        s = '%3d digits, %-40s: %d primes, proved %d, no chain %d, timed out (%d s) %d, checker rejects %d, median steps %s, median time %s s, longest %s s' % (
            d, label, count, ok, fail, tlim, tout, bad, row['median_steps'], '%.2f' % row['median_s'] if times else '-', '%.2f' % row['max_s'] if times else '-')
        print(s, flush=True)
        lines.append(s)
    if pari_ok is not None:
        s = '%3d digits: PARI isprime (a proof) agrees on %d of %d' % (d, pari_ok, count)
        print(s, flush=True)
        lines.append(s)
json.dump(rows, open('chain_test_%s.json' % tag, 'w'))
open('chain_test_%s.txt' % tag, 'w').write('\n'.join(lines) + '\n')
