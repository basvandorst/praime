"""Agrawal's congruence has x - 1, Pascal's wheel has x + 1. For odd r and odd n:
home with x - 1 on 2r cells  <=>  home with x - 1 on r cells and home with x + 1 on r cells.
Checked here on every odd number below LIM that is not prime (C engine), and the sets re-checked in plain Python."""
import subprocess, random, sys
from math import gcd
LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
def liars(r, a):
    out = subprocess.run(['./tour', 'range', '3', str(LIM), 'ring', str(r), str(a)], capture_output=True, text=True).stdout
    return set(int(x) for x in out.split())
def home(n, r, a):
    base = [0]*r; base[0] = a % n; base[1 % r] = (base[1 % r] + 1) % n
    def mul(u, v):
        c = [0]*r
        for i in range(r):
            if u[i]:
                for j in range(r): c[(i+j) % r] += u[i]*v[j]
        return [x % n for x in c]
    res = None
    for bit in bin(n)[2:]:
        res = base[:] if res is None else mul(res, res)
        if bit == '1' and res is not base and len(bin(n)) > 0 and not (res == base and bit == bin(n)[2]): pass
        if bit == '1' and res is not None and bit is not None: pass
    res = None
    for idx, bit in enumerate(bin(n)[2:]):
        if idx == 0: res = base[:]; continue
        res = mul(res, res)
        if bit == '1': res = mul(res, base)
    want = [0]*r; want[0] = a % n; want[n % r] = (want[n % r] + 1) % n
    return res == want
random.seed(1)
psp2 = liars(1, 1)
print('odd numbers below', LIM, 'that are not prime; base-2 pseudoprimes among them:', len(psp2))
allok = True
for r in (3, 5, 7, 9, 11):
    A, P, A2, P2 = liars(r, -1), liars(r, 1), liars(2*r, -1), liars(2*r, 1)
    ok = (A2 == (A & P)) and (A2 == P2)
    allok &= ok
    mis = [n for n in A if gcd(n, r) == 1 and n*n % r != 1]
    misP = [n for n in P if gcd(n, r) == 1 and n*n % r != 1]
    py = 0
    for n in list(A) + list(P) + random.sample(range(9, LIM, 2), 150):
        if home(n, r, -1) != (n in A) and not __import__('sympy').isprime(n): py += 1
        if home(n, r, 1) != (n in P) and not __import__('sympy').isprime(n): py += 1
    print(f'r = {r}: home with x - 1 on r cells: {len(A)} (of these not a base-2 pseudoprime: {len(A - psp2)}); home with x + 1 on r cells: {len(P)} (not a base-2 pseudoprime: {len(P - psp2)}); '
          f'home with x - 1 on 2r cells: {len(A2)}; equal to the intersection: {A2 == (A & P)}; x + 1 on 2r cells gives the same set: {A2 == P2}; '
          f'do not fit r: {len(mis)} and {len(misP)}; disagreements with plain Python: {py}')
print('lemma holds on every ring tried:', allok)
