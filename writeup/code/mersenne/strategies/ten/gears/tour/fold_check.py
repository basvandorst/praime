"""The folded wheel against plain Python, and the proof-by-data that folding loses no liar.
Fold of the ring of r cells (r prime) by the subgroup H of e-th powers mod r (e divides (r-1)/2):
u = product over h in H of (x^h + 1); home on r cells implies u^n = u with x -> x^n (mod n, x^r - 1)."""
import json, random, subprocess, sys
from math import gcd
def prim_root(r):
    for g in range(2, r):
        x, ok = 1, True
        for k in range(1, r-1):
            x = x*g % r
            if x == 1: ok = False; break
        if ok: return g
def classes(r, e):
    g = prim_root(r); ind = {}; x = 1
    for k in range(r-1): ind[x] = k; x = x*g % r
    return {x: ind[x] % e for x in range(1, r)}
def conv(a, b, r, n):
    c = [0]*r
    for i in range(r):
        if a[i]:
            for j in range(r):
                if b[j]: c[(i+j) % r] += a[i]*b[j]
    return [v % n for v in c]
def fold_vec(n, r, e):
    cl = classes(r, e); P = [0]*r; P[0] = 1
    for h in range(1, r):
        if cl[h] == 0:
            f = [0]*r; f[0] = 1; f[h] += 1; P = conv(P, f, r, n)
    res = None
    for i, bit in enumerate(bin(n)[2:]):
        if i == 0: res = P[:]; continue
        res = conv(res, res, r, n)
        if bit == '1': res = conv(res, P, r, n)
    want = [0]*r
    for z in range(r): want[z*n % r] = P[z]
    rep = {}
    for x in range(r-1, 0, -1): rep[cl[x]] = x
    cells = [(res[rep[k]] - res[0]) % n for k in range(e)]
    ok = all((res[z] - res[0] - want[z] + want[0]) % n == 0 for z in range(r))
    return ok, cells
def folds_of(r):
    m = (r-1)//2; out = []
    for q in range(2, m+1):
        if m % q == 0 and all(q % d for d in range(2, q)):
            e = q
            while m % e == 0: out.append(e); e *= q
    return sorted(out)
if __name__ == "__main__":
    random.seed(20261118)
    primes = [5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 43]
    psp = [int(x) for x in open('psp2_to1e10.txt')]
    big = [int(x) for x in open('sample64.txt')]
    txt = '\n'.join(map(str, psp)) + '\n'
    lo = subprocess.run(['./ring3', '-plain', '-fit'] + [str(r) for r in primes], input=txt, capture_output=True, text=True).stdout
    liars = {}
    for l in lo.split('\n'):
        if l.startswith('ring'): _, r, _, n = l.split(); liars.setdefault(int(r), set()).add(int(n))
    al = subprocess.run(['./fold', '-all'] + [str(r) for r in primes], input=txt, capture_output=True, text=True).stdout
    passed = {}; nfold = 0
    for l in al.split('\n'):
        f = l.split()
        if len(f) == 4:
            nfold += 1
            if f[0] == 'P': passed.setdefault((int(f[1]), int(f[2])), set()).add(int(f[3]))
    lost = 0; tot = 0
    for r in primes:
        for e in folds_of(r):
            for n in liars.get(r, ()):
                if True:
                    tot += 1
                    if n not in passed.get((r, e), ()): lost += 1
    print('rings', primes, '| folds per ring:', {r: folds_of(r) for r in primes})
    print('liars on the full ring below 1e10 (numbers that fit), times folds:', tot, '| liars that a fold would have thrown away:', lost)
    for r in primes:
        print('  ring %d: %d liars on the full ring;' % (r, len(liars.get(r, ()))), ', '.join('fold of %d cells lets %d of %d candidates through' % (e, len(passed.get((r, e), ())), sum(1 for n in psp if gcd(n, r) == 1)) for e in folds_of(r)))
    pick = random.sample(big, 25) + [n for n in big if n > 2**63][:10] + random.sample(psp, 25) + [18446744073709551557, 18446744073709551615]
    for r in primes: pick += random.sample(sorted(liars.get(r, {561})), min(3, len(liars.get(r, {561}))))
    du = subprocess.run(['./fold', '-all', '-dump'] + [str(r) for r in primes], input='\n'.join(map(str, pick)) + '\n', capture_output=True, text=True).stdout.split('\n')
    cmpd = wrong = 0; verdict = {}
    for l in du:
        f = l.split()
        if len(f) == 4: verdict[(int(f[3]), int(f[1]), int(f[2]))] = f[0] == 'P'
    for l in du:
        f = l.split()
        if f and f[0] == 'D':
            n, r, e = int(f[1]), int(f[2]), int(f[3]); ok, cells = fold_vec(n, r, e); cmpd += 1
            if cells != [int(x) for x in f[4:]] or ok != verdict[(n, r, e)]: wrong += 1
    print('folded wheels compared with plain Python on the full ring (cells and verdict), numbers up to 2^64:', cmpd, '| different:', wrong)
    primes2 = primes + [23]
    o = subprocess.run(['./fold', '-fit'] + [str(r) for r in primes2], input=txt, capture_output=True, text=True)
    got = {}
    for l in o.stdout.split('\n'):
        if l.startswith('ring'): _, r, _, n = l.split(); got.setdefault(int(r), set()).add(int(n))
    lo23 = subprocess.run(['./ring3', '-plain', '-fit', '23'], input=txt, capture_output=True, text=True).stdout
    liars[23] = set(int(l.split()[3]) for l in lo23.split('\n') if l.startswith('ring'))
    print('whole pipeline (fold or half ring, then the full ring) on every candidate below 1e10, numbers that fit included: same liars as the plain full ring on', sum(1 for r in primes2 if got.get(r, set()) == liars.get(r, set())), 'of', len(primes2), 'rings;', sum(len(v) for v in got.values()), 'liars')
    for name, lst in (('candidates below 1e10, numbers that do not fit', psp), ('sample of the list up to 2^64, numbers that do not fit', big)):
        o = subprocess.run(['./fold'] + [str(r) for r in primes2], input='\n'.join(map(str, lst)) + '\n', capture_output=True, text=True)
        print(name + ':'); print(o.stderr.strip()); assert not o.stdout.strip()
