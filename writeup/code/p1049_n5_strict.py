"""Puzzle 1049, n = 5, strict reading (all chain members zero-free primes).
p -> q = f(p) must be a zero-free prime with chain(q) >= 4.  Since q <= 9^9 * digits(p) and the least
such q is the 12-digit 112222268681, every candidate q is an f-image of a 12-digit multiset, and
p1049_search.py showed that 111222226688 is the only multiset of <= 12 digits whose image has chain >= 3:
so q runs over the prime arrangements of 111222226688 (5566 of them).  For each q, D(q) = least
number of values d^d (d = 1..9) summing to q (coin change); then the least p is the smallest prime
arrangement of an optimal multiset, over all q with the least D.   Output: data/p1049_n5_strict.txt"""
import time
from gmpy2 import mpz, is_prime
from sympy.utilities.iterables import multiset_permutations
POW = [d**d for d in range(10)]
def f(x): return sum(POW[int(c)] for c in str(x))
def chain(p):
    k = 0
    while k < 20:
        q = f(p)
        if not is_prime(q) or '0' in str(q): return k
        k += 1; p = q
    return k
def min_digits(q):
    best = [None]; sols = []
    def rec(d, rem, cnt, used):
        if d == 0:
            if rem == 0:
                if best[0] is None or cnt < best[0]: best[0] = cnt; sols.clear(); sols.append(tuple(used))
                elif cnt == best[0]: sols.append(tuple(used))
            return
        if best[0] is not None and cnt >= best[0]: return
        # lower bound: remaining value needs at least ceil(rem / POW[d]) more digits
        if best[0] is not None and cnt + -(-rem // POW[d]) >= best[0]: return
        top = rem // POW[d]
        for c in range(top, -1, -1):
            used.append(c); rec(d - 1, rem - c * POW[d], cnt + c, used); used.pop()
    rec(9, q, 0, [])
    return best[0], sols
t0 = time.time()
qs = sorted(int("".join(p)) for p in multiset_permutations(list("111222226688")))
qs = [q for q in qs if is_prime(q)]
assert len(qs) == 5566 and all(chain(q) >= 4 for q in qs)
res = []
for q in qs:
    D, sols = min_digits(q); res.append((D, q, sols))
res.sort()
Dmin = res[0][0]
hist = {}
for D, q, s in res: hist[D] = hist.get(D, 0) + 1
print("D(q) histogram over the 5566 candidates:", dict(sorted(hist.items())), f"{time.time()-t0:.0f}s")
win = [(q, sols) for D, q, sols in res if D == Dmin]
print(f"least digit count {Dmin}, attained by {len(win)} q:", [q for q, _ in win])
out = open("writeup/data/p1049_n5_strict.txt", "w")
out.write(f"# puzzle 1049 n=5 strict: D(q) = least digit count of a zero-free p with f(p) = q, q over the prime arrangements of 111222226688\n")
for D, q, sols in res: out.write(f"q={q} D={D} multisets={['9^%d 8^%d 7^%d 6^%d 5^%d 4^%d 3^%d 2^%d 1^%d' % s for s in sols]}\n")
bestp = None
for q, sols in win:
    for counts in sols:
        digits = []
        for d in range(1, 10): digits += [str(d)] * counts[9 - d]
        arr = sorted(digits); n = 0
        while True:
            n += 1
            if bestp is not None and "".join(arr) >= bestp: break
            if is_prime(mpz("".join(arr))):
                bestp = "".join(arr); out.write(f"candidate p={bestp} (q={q}, arrangement number {n})\n")
                print(f"q={q}: smallest prime arrangement found after {n} arrangements: {bestp[:40]}...{bestp[-30:]}"); break
            i = len(arr) - 2
            while i >= 0 and arr[i] >= arr[i + 1]: i -= 1
            if i < 0: break
            j = len(arr) - 1
            while arr[j] <= arr[i]: j -= 1
            arr[i], arr[j] = arr[j], arr[i]; arr[i + 1:] = reversed(arr[i + 1:])
print(f"least p for n = 5 (strict): {len(bestp)} digits\n{bestp}")
out.write(f"least p: {bestp}\n"); out.close()
