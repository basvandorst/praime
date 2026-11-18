import json, random, sys, time
from chain_wheel import certify, sprp, SMALL
from chain_verify import verify
import sympy

out = {}
lines = []


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
for label, kw in (('every size of wheel', {}), ('curves only', {'use_pm': False})):
    t = time.time()
    proved = wrong = undecided = bad = comp_into = 0
    kinds = {}
    und = []
    for n in range(2, TOP + 1):
        cert, nodes = certify(n, stop=50, **kw)
        isp = sympy.isprime(n)
        if isinstance(cert, list):
            if not isp:
                wrong += 1
            if not verify(n, cert):
                bad += 1
            proved += 1
            for s in cert:
                kinds[s['t']] = kinds.get(s['t'], 0) + 1
        else:
            if isp:
                undecided += 1
                und.append(n)
            elif sprp(n):
                comp_into += 1
    say('chain with the proof stopped only below 50,', label, ': numbers 2 to', TOP, ': primes proved', proved, ', primes without a chain', undecided, und[:20], ', composites given a chain', wrong, ', chains the checker rejects', bad, ', composites that passed the doubling wheel and entered the search', comp_into, ', steps by kind', sorted(kinds.items()), '(%.0f s)' % (time.time() - t))
    out[label] = dict(top=TOP, proved=proved, undecided=undecided, und=und, wrong=wrong, bad=bad, comp_into=comp_into, kinds=kinds)
json.dump(out, open('chain_test_small.json', 'w'))
open('chain_test_small.txt', 'w').write('\n'.join(lines) + '\n')
