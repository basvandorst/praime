import json, sys
from chain_verify import verify
import sympy

rows = json.load(open(sys.argv[1]))
label = sys.argv[2]
known = sys.argv[3] if len(sys.argv) > 3 else 'check'


def conv(st):
    o = {}
    for k, v in st.items():
        if k == 't':
            o[k] = v
        elif k == 'f':
            o[k] = [[int(f), int(a)] for f, a in v]
        else:
            o[k] = int(v)
    return o


proved = und = wrong = bad = comp = missed = 0
kinds = {}
ms = []
for row in rows:
    n = int(row[0])
    ok, u, steps = row[1], row[2], row[3]
    isp = sympy.isprime(n) if known == 'check' else known == 'prime'
    if ok:
        proved += 1
        cert = [conv(s) for s in steps]
        if not isp:
            wrong += 1
        if not verify(n, cert):
            bad += 1
        for s in cert:
            kinds[s['t']] = kinds.get(s['t'], 0) + 1
    elif isp:
        und += 1
    else:
        comp += 1
    if len(row) > 4:
        ms.append(row[4])
s = '%s: %d numbers, proved prime %d, primes without a chain %d, composites given a chain %d, chains the Python checker rejects %d, composites turned away %d, steps by kind %s' % (
    label, len(rows), proved, und, wrong, bad, comp, sorted(kinds.items()))
if ms:
    ms.sort()
    s += ', median time %d ms, longest %d ms' % (ms[len(ms) // 2], ms[-1])
print(s)
open(sys.argv[4], 'a').write(s + '\n') if len(sys.argv) > 4 else None
