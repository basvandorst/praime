"""Independent recheck of the derived real sequences (plain integers, sympy, own loops; nothing from gen.py).
  pi(a_n): sympy.primepi.  Class counts: a plain (non-segmented) odd-only sieve to B with its own code, checked
  against pi.  TF survivors < 10^6: Python pow(2, p, q) for every q = 2jp+1, j = 1..1000 (no mod-8 shortcut),
  q < 2^p - 1.  Bit/decimal-length indices: plain int arithmetic + sympy.primepi.
  Null lists: sorted, distinct, prime (sympy.isprime on all), 51 terms; inclusion frequencies of small primes.
"""
import numpy as np, sympy, math, time
D = np.load('data.npz')
REAL = [int(x) for x in D['real']]
B = 83_880_077
ok = True
t0 = time.time()
pi = [int(sympy.primepi(a)) for a in REAL]
print('pi via sympy.primepi: %.1f s' % (time.time() - t0))
s = np.ones(B // 2, bool); s[0] = False      # s[i] <-> 2i+1
for i in range(1, (math.isqrt(B) + 1) // 2 + 1):
    if s[i]:
        p = 2 * i + 1; s[p * p // 2::p] = False
odd = 2 * np.nonzero(s)[0] + 1
del s
print('pi(B) =', len(odd) + 1, '(sympy:', int(sympy.primepi(B)), ')')
n1 = np.cumsum(odd % 4 == 1); n3 = np.cumsum(odd % 4 == 3)


def cls_index(a):
    k = int(np.searchsorted(odd, a))
    assert odd[k] == a
    return int(n1[k]) if a % 4 == 1 else int(n3[k])


der = {}
der['pi'] = pi
der['gap'] = [pi[i + 1] - pi[i] for i in range(50)]
der['idx1'] = [cls_index(a) for a in REAL if a % 4 == 1]
der['idx3'] = [cls_index(a) for a in REAL if a % 4 == 3]
der['cls'] = [cls_index(a) for a in REAL if a > 2]
der['blen'] = [a - 2 ** (a.bit_length() - 1) + 1 for a in REAL]
der['pblen'] = [int(sympy.primepi(a)) - int(sympy.primepi(2 ** (a.bit_length() - 1) - 1)) for a in REAL]
der['pdec'] = [int(sympy.primepi(a)) - int(sympy.primepi(10 ** (len(str(a)) - 1) - 1)) for a in REAL]
# TF survivors
t0 = time.time()
surv = []
for p in sympy.primerange(2, 10 ** 6):
    M = (1 << p) - 1 if p < 64 else None
    has = False
    for j in range(1, 1001):
        q = 2 * j * p + 1
        if M is not None and q >= M: break
        if pow(2, p, q) == 1:
            has = True; break
    if not has: surv.append(p)
print('TF survivors (plain loop): %d, %.1f s' % (len(surv), time.time() - t0))
ss = {v: i + 1 for i, v in enumerate(surv)}
open('check_surv.txt', 'w').write(' '.join(map(str, surv)) + '\n')
der['tf'] = [ss[a] for a in REAL if a < 10 ** 6]
assert len(surv) == len(D['surv']) and all(int(a) == b for a, b in zip(D['surv'], surv)), 'survivor list differs'
for nm in [str(x) for x in D['names']]:
    mine = [int(x) for x in D[nm][0, :D[nm + '_len'][0]]]
    same = mine == der[nm]
    ok &= same
    print(f'{nm:6s} {"OK " if same else "DIFF"} {len(mine)} terms; first 15: {der[nm][:15]}')
# null lists
NW = D['nullW']; NT = D['nullT']
allp = set()
for row in NW:
    r = [int(x) for x in row]
    assert len(r) == 51 and r == sorted(set(r)) and r[-1] < B
    allp.update(r)
for row in NT:
    r = [int(x) for x in row]
    assert len(r) == 33 and r == sorted(set(r)) and r[-1] < 10 ** 6 and all(x in ss for x in r)
assert all(sympy.isprime(x) for x in allp)
print('null W: %d lists, %d distinct primes, all prime; null TF: %d lists of survivors' % (len(NW), len(allp), len(NT)))
EG = math.exp(0.5772156649015329)
for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 61, 89):
    a = 6 if p % 4 == 1 else 2
    print(f'   P({p:2d} in list): model {min(1, EG * math.log2(a * p) / p):.3f}, empirical {np.mean((NW == p).any(1)):.3f}')
print('ALL DERIVED SEQUENCES REPRODUCED' if ok else 'MISMATCH')
