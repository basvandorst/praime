# Puzzle 360: N such that for every divisor d, d + N/d and |d - N/d| are 1 or prime.  (about 30 s with L = 5000)
# Claims (results/360.md): the valid N <= 10^5 are 2, 6, 30, 42, 462; every valid N > 462 is divisible by 2, 3, 7, 11, 19
# (cascade from {2, 3}, the only exceptions being 30, 42, 462) and by every prime = 3 (mod 4) below L.
# Forcing rule: for a prime l not dividing N and d | N, l | d - N/d iff N = d^2 (mod l), l | d + N/d iff N = -d^2 (mod l);
# a value divisible by l must equal l; d + N/d >= 2 sqrt(N) > l (N > 10^5, l < 632, or N >= the forced product), and
# |d - N/d| = l for at most one pair. So if every r != 0 (mod l) is +-d^2 for >= 3 distinct subset products d of the
# forced primes, l | N; residues reached by fewer d leave the candidates N = d(d +- l), which are tested.
# Run: .venv/bin/python results/code/p360_cascade.py 5000 --small
import sys, time
from sympy import primerange, isprime, divisors
def valid(N): return all((x == 1 or isprime(x)) for d in divisors(N) for x in (d + N // d, abs(d - N // d)))
# (i) every valid N <= 10^5, directly
small = [N for N in range(1, 10**5 + 1) if valid(N)] if '--small' in sys.argv else [2, 6, 30, 42, 462]
print("valid N <= 1e5:", small)
def test(F, l):
    # cnt[x] = number of distinct subset products d of F with d = x (mod l), capped at 3; rep[x] = those d (when < 3)
    cnt = [0] * l; rep = [[] for _ in range(l)]; cnt[1] = 1; rep[1] = [1]
    for p in F:
        pm = p % l
        ncnt = cnt[:]; nrep = [r[:] for r in rep] if True else None
        for x in range(1, l):
            if cnt[x]:
                y = x * pm % l
                if ncnt[y] < 3:
                    for d in rep[x]:
                        if ncnt[y] < 3: nrep[y].append(d * p); ncnt[y] += 1
                    ncnt[y] = min(3, ncnt[y] + max(0, cnt[x] - len(rep[x])))
        cnt, rep = ncnt, nrep
    # residue r of N is excluded by every d with d^2 = +-r; l = 3 (mod 4): each r != 0 is +x^2 or -x^2 for x = +-sqrt
    need = []; hits = [0] * l; hd = [[] for _ in range(l)]
    for x in range(1, l):
        if cnt[x]:
            for r in {x * x % l, (-x * x) % l}:
                hits[r] += cnt[x]; hd[r] += rep[x]
    for r in range(1, l):
        if hits[r] == 0: return False, None
        if hits[r] < 3: need += hd[r]
    return True, need
F = [2, 3]; prod = 6                       # N > 6 valid: N - 1, N + 1 prime => 6 | N
for q, exc in ((7, [30]), (11, [42]), (19, [462])):
    ok, need = test(F, q); assert ok
    cands = sorted({N for d in need for N in (d * (d + q), d * (d - q)) if N > 6 and N % prod == 0})
    assert [N for N in cands if valid(N)] == exc, (q, cands)    # the only escape from q | N is a known term
    F.append(q); prod *= q
assert small == [2, 6, 30, 42, 462]
L = int(sys.argv[1]); pend = []; t0 = time.time()
for l in primerange(23, L):
    if l % 4 != 3: continue
    pend.append(l)
    progress = True
    while progress:
        progress = False
        for q in list(pend):
            ok, need = test(F, q)
            if not ok: continue
            cands = [N for d in need for N in (d * (d + q), d * (d - q)) if N > 10**5 and N % prod == 0]
            assert not any(valid(N) for N in cands)
            F.append(q); prod *= q; pend.remove(q); progress = True
print("not forced (pending):", pend[:20], len(pend))
print("forced:", len(F), "primes; last", F[-5:], "; product has", len(str(prod)), "digits; time", round(time.time() - t0), "s")
allq = [q for q in primerange(3, L) if q % 4 == 3]
print("all primes = 3 (mod 4) below", L, "forced:", set(allq) <= set(F))
