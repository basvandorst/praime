"""Puzzle 282: primes p with phi(p) = phi(rev p) or sigma(p) = sigma(rev p).

Write r = rev(p) = m*b with b prime, b not dividing m. Then
  phi:   m*rev(r) = phi(m)*r   - m*(phi(m) - 1)
  sigma: m*rev(r) = sigma(m)*r + m*(sigma(m) - 1)
a linear "reverse multiple" equation P*rev(r) = Q*r + C. A two-ended carry
automaton (low carry from the units digit, high carry from the leading digit)
finds every solution r of every length; then b = r/m and p = rev(r) are tested.

Usage: p282_search.py phi|sigma M KMIN KMAX
Complete for: every m <= M with 3, 11, 10 not dividing m, every length k in
[KMIN, KMAX] with 10^ceil(k/2) > |C| (shorter k are skipped and reported).
"""
import sys
import gmpy2
from sympy import totient, divisor_sigma


def levels_of(P, Q, C, kmax):
    B = max(P, Q) + 1
    levels = [{(-C, 0): []}]
    for t in range(kmax // 2 + 1):
        nxt = {}
        xs = range(1, 10) if t == 0 else range(10)
        ys = (1, 3, 7, 9) if t == 0 else range(10)
        for cl, ch in levels[-1]:
            for x in xs:
                for y in ys:
                    s = cl + P * y - Q * x
                    if s % 10:
                        continue
                    nl = s // 10
                    if abs(nl) > B + abs(C) // 10 ** (t + 1) + 1:
                        continue
                    nh = 10 * ch - (P * x - Q * y)
                    if abs(nh) > B:
                        continue
                    nxt.setdefault((nl, nh), []).append(((cl, ch), x, y))
        levels.append(nxt)
    return levels


def paths(levels, t, state):
    if t == 0:
        yield [], []
        return
    for prev, x, y in levels[t][state]:
        for lo, hi in paths(levels, t - 1, prev):
            yield lo + [x], hi + [y]


def solutions(P, Q, C, k, levels):
    t = k // 2
    for (cl, ch) in levels[t]:
        if k % 2:
            mids = [d for d in range(10) if cl + (P - Q) * d == 10 * ch]
        else:
            mids = [None] if cl == ch else []
        for d in mids:
            for lo, hi in paths(levels, t, (cl, ch)):
                digs = lo + ([] if d is None else [d]) + hi[::-1]
                yield int(''.join(map(str, reversed(digs))))


def main():
    kind, M, kmin, kmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    skipped = 0
    for m in range(2, M + 1):
        if m % 3 == 0 or m % 11 == 0 or m % 10 == 0:
            continue
        f = int(totient(m)) if kind == 'phi' else int(divisor_sigma(m))
        C = -m * (f - 1) if kind == 'phi' else m * (f - 1)
        levels = levels_of(m, f, C, kmax)
        for k in range(kmin, kmax + 1):
            if 10 ** ((k + 1) // 2) <= abs(C):
                skipped += 1
                continue
            for r in solutions(m, f, C, k, levels):
                if r % m:
                    continue
                b = r // m
                if m % b == 0 or not gmpy2.is_prime(b):
                    continue
                p = int(str(r)[::-1])
                if p != r and gmpy2.is_prime(p):
                    print(kind, 'm=%d' % m, 'digits=%d' % k, 'p=%d' % p, flush=True)
    print('# done', kind, 'M =', M, 'k =', kmin, '..', kmax, '; (m, k) pairs skipped for small k:', skipped)


main()
