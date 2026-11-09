"""Puzzle 1049 under the loose reading: only p must be zero-free; later chain members may contain zeros.
n = 4: multisets of <= 12 digits whose image v is prime with loose chain >= 3; least prime arrangement.
n = 5: q = f(p) must be a prime with loose chain >= 4 and q <= 9^9 D for the digit count D of p; the
least such q has 12 digits (checked below), and q < 400 * 9^9 < 1.6e11, so q is a prime arrangement of
an 11- or 12-digit multiset (zeros allowed in q) whose image has loose chain >= 3.
D(q) = least number of values d^d summing to q; least p = smallest prime arrangement.
Output: data/p1049_loose.txt"""
import itertools, time
from gmpy2 import mpz, is_prime
from sympy.utilities.iterables import multiset_permutations
POW = [d**d for d in range(10)]
def f(x): return sum(POW[int(c)] for c in str(x))
def lchain(p):
    k = 0
    while k < 20:
        q = f(p)
        if not is_prime(q): return k
        k += 1; p = q
    return k
t0 = time.time()
out = open("writeup/data/p1049_loose.txt", "w")
qual = []
for size in range(1, 13):
    for ms in itertools.combinations_with_replacement(range(1, 10), size):
        v = sum(POW[d] for d in ms)
        if is_prime(v) and lchain(v) >= 3: qual.append(ms)
print("n=4 loose: qualifying zero-free multisets (<= 12 digits):", ["".join(map(str, m)) for m in qual], flush=True)
best4 = None
for ms in qual:
    for perm in multiset_permutations(list(ms)):
        v = int("".join(map(str, perm)))
        if is_prime(v): best4 = v if best4 is None else min(best4, v); break
print("n=4 loose least p:", best4, "chain:", lchain(best4)); out.write(f"n=4 loose least p: {best4}\n")
cands = []
for size in (11, 12):
    for ms in itertools.combinations_with_replacement(range(0, 10), size):
        if ms[-1] == 0: continue
        v = sum(POW[d] for d in ms if d)
        if is_prime(v) and lchain(v) >= 3:
            for perm in multiset_permutations(list(ms)):
                if perm[0] == 0: continue
                q = int("".join(map(str, perm)))
                if is_prime(q) and q < 400 * POW[9]: cands.append(q)
cands = sorted(set(cands)); print(len(cands), "candidate q (11-12 digits, loose chain >= 4, q < 400 * 9^9);", f"{time.time()-t0:.0f}s", flush=True)
assert all(lchain(q) >= 4 for q in cands)
def min_digits(q):
    best = [None]; sols = []
    def rec(d, rem, cnt, used):
        if d == 0:
            if rem == 0:
                if best[0] is None or cnt < best[0]: best[0] = cnt; sols.clear(); sols.append(tuple(used))
                elif cnt == best[0]: sols.append(tuple(used))
            return
        if best[0] is not None and cnt + -(-rem // POW[d]) >= best[0]: return
        for c in range(rem // POW[d], -1, -1):
            used.append(c); rec(d - 1, rem - c * POW[d], cnt + c, used); used.pop()
    rec(9, q, 0, [])
    return best[0], sols
res = sorted((min_digits(q)[0], q) for q in cands)
Dmin = res[0][0]; win = [q for D, q in res if D == Dmin]
print("least digit count", Dmin, "for q in", win, flush=True)
out.write(f"# n=5 loose: D(q) for the {len(cands)} candidates\n" + "".join(f"q={q} D={D}\n" for D, q in res))
bestp = None
for q in win:
    for counts in min_digits(q)[1]:
        digits = []
        for d in range(1, 10): digits += [str(d)] * counts[9 - d]
        arr = sorted(digits); n = 0
        while True:
            n += 1
            if bestp is not None and "".join(arr) >= bestp: break
            if is_prime(mpz("".join(arr))): bestp = "".join(arr); out.write(f"candidate p={bestp} (q={q}, arrangement {n})\n"); break
            i = len(arr) - 2
            while i >= 0 and arr[i] >= arr[i + 1]: i -= 1
            if i < 0: break
            j = len(arr) - 1
            while arr[j] <= arr[i]: j -= 1
            arr[i], arr[j] = arr[j], arr[i]; arr[i + 1:] = reversed(arr[i + 1:])
print(f"n=5 loose least p: {len(bestp)} digits, f(p) = {f(int(bestp))}, chain {lchain(int(bestp))}")
print(bestp); out.write(f"least p: {bestp}\n"); out.close()
