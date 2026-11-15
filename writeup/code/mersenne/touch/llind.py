"""Independent Lucas-Lehmer check (written separately from ../ll.py): reduction mod M = 2^p - 1 by shift and add,
primes from a plain Python sieve (no sympy).  Modes:
  llind.py all X        -> LL on every prime p < X, print the exponents found
  llind.py list p1,p2.. -> LL on the given p, print p and the verdict
  llind.py time p1,p2.. -> wall time of one LL per p (cost model)"""
import sys, time, gmpy2

def primes_below(n):
    s = bytearray([1]) * n; s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(range(i * i, n, i)))
    return [i for i in range(n) if s[i]]

def ll(p):
    if p == 2: return True
    M = (gmpy2.mpz(1) << p) - 1; s = gmpy2.mpz(4)
    for _ in range(p - 2):
        s = s * s - 2
        s = (s & M) + (s >> p)
        if s >= M: s -= M
    return s == 0

if __name__ == '__main__':
    mode, arg = sys.argv[1], sys.argv[2]
    if mode == 'all':
        X = int(arg); t = time.time()
        out = [p for p in primes_below(X) if ll(p)]
        print(f'LL on all primes < {X}: {len(out)} exponents {out}  ({time.time() - t:.0f} s)')
    else:
        for p in map(int, arg.split(',')):
            t = time.time(); r = ll(p); dt = time.time() - t
            print(p, 'PRIME' if r else 'composite', f'{dt:.3f} s', f'{dt / max(p - 2, 1) * 1e6:.2f} us/squaring')
