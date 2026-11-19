"""The gaps between consecutive Mersenne exponents, tested against two random models.

Data: the 52 known exponents. GIMPS has verified every exponent below 83,200,483 (FACTS.md, fact 6301), so the
first 51 are certainly consecutive; the gap from 82,589,933 to 136,279,841 is reported apart ("52?").

Model for the sizes (Lenstra, Pomerance, Wagstaff; a conjecture): log2 of the exponents behaves like a Poisson
process with rate e^gamma per doubling, so log2(p_(n+1)/p_n) is exponential with mean e^(-gamma) = 0.5615.

Models for the remainders, used for "is the gap a multiple of m" (equivalently: p_n and p_(n+1) lie in the same
class mod m). Each real exponent p >= 61 is replaced by a random prime of about the same size:
  A: any prime;
  B: a prime for which 2^q - 1 has no factor 2kq + 1 with k <= 100 (a Mersenne exponent has none at all).
The window is [p/1.5, 1.5p] below 10,000 and [0.95p, 1.05p] above.

usage: mexp_gaps.py SIMS      writes mexp_gaps.json and mexp_gaps.txt"""
import sys, json, math, random
import gmpy2

E = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213,
     19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377,
     6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161,
     74207281, 77232917, 82589933, 136279841]
CERTAIN = 51
GAMMA = 0.5772156649015329
K = 100
MODS = [3, 4, 5, 7, 8, 9, 11, 12, 13, 16, 24]
WHEELS = list(range(2, 61))


def survives(q, K=K):
    for k in range(1, K + 1):
        f = 2 * k * q + 1
        if f % 8 in (1, 7) and gmpy2.is_prime(f) and pow(2, q, f) == 1:
            return False
    return True


def rand_prime(rng, p, weighted):
    lo, hi = (int(p / 1.5), int(p * 1.5)) if p < 10000 else (int(p * 0.95), int(p * 1.05))
    lo = max(lo, 3)
    while True:
        q = rng.randrange(lo, hi + 1) | 1
        if gmpy2.is_prime(q) and (not weighted or survives(q)):
            return q


def gaps_divisible(seq, m):
    return sum(1 for a, b in zip(seq, seq[1:]) if (b - a) % m == 0)


def main():
    sims = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    out = {}
    lines = []

    def say(s):
        lines.append(s)
        print(s)
        sys.stdout.flush()

    assert len(E) == 52 and all(gmpy2.is_prime(p) for p in E) and E == sorted(E)
    ex = E[:CERTAIN]
    g = [b - a for a, b in zip(ex, ex[1:])]
    out['exponents'] = E
    out['gaps'] = [b - a for a, b in zip(E, E[1:])]
    say('gaps between the first 51 exponents (certainly consecutive): ' + ', '.join(str(x) for x in g))
    say('gap to the 52nd (not yet verified to be the next): %d' % (E[51] - E[50]))
    say('even gaps: %d of %d (only 2 -> 3 is odd)' % (sum(1 for x in g if x % 2 == 0), len(g)))

    # 1. sizes: log2 ratios against the exponential law with mean e^-gamma
    lam = math.exp(GAMMA)
    for start in (5, 127):
        seq = [p for p in ex if p >= start]
        r = [math.log2(b / a) for a, b in zip(seq, seq[1:])]
        n = len(r)
        mean = sum(r) / n
        u = sorted(1 - math.exp(-lam * x) for x in r)
        D = max(max((i + 1) / n - u[i], u[i] - i / n) for i in range(n))
        rho = sum((r[i] - mean) * (r[i + 1] - mean) for i in range(n - 1)) / sum((x - mean) ** 2 for x in r)
        rng = random.Random(1)
        cD = cM = cR = 0
        for _ in range(20000):
            z = [rng.expovariate(lam) for _ in range(n)]
            mz = sum(z) / n
            uz = sorted(1 - math.exp(-lam * x) for x in z)
            Dz = max(max((i + 1) / n - uz[i], uz[i] - i / n) for i in range(n))
            rz = sum((z[i] - mz) * (z[i + 1] - mz) for i in range(n - 1)) / sum((x - mz) ** 2 for x in z)
            cD += Dz >= D
            cM += abs(mz - 1 / lam) >= abs(mean - 1 / lam)
            cR += abs(rz) >= abs(rho)
        out['sizes_from_%d' % start] = dict(n=n, mean=mean, model=1 / lam, D=D, pD=cD / 20000, pMean=cM / 20000, rho=rho, pRho=cR / 20000,
                                            largest=max(r), smallest=min(r))
        say('sizes, exponents from %d on (%d gaps): mean log2 ratio %.3f (model %.3f, p = %.2f); KS distance %.3f (p = %.2f); lag-1 correlation %.2f (p = %.2f); largest ratio 2^%.2f, smallest 2^%.4f'
            % (start, n, mean, 1 / lam, cM / 20000, D, cD / 20000, rho, cR / 20000, max(r), min(r)))
    # counts per doubling: dispersion
    lo, hi = 2, 26
    cnt = [sum(1 for p in ex if 2 ** k <= p < 2 ** (k + 1)) for k in range(lo, hi)]
    mu = sum(cnt) / len(cnt)
    var = sum((c - mu) ** 2 for c in cnt) / (len(cnt) - 1)
    rng = random.Random(2)
    cV = 0
    for _ in range(20000):
        z = []
        for k in range(lo, hi):
            t, c = 0.0, 0
            while True:
                t += rng.expovariate(mu)
                if t > 1:
                    break
                c += 1
            z.append(c)
        mz = sum(z) / len(z)
        vz = sum((c - mz) ** 2 for c in z) / (len(z) - 1)
        cV += vz / max(mz, 1e-9) <= var / mu
    out['doublings'] = dict(lo=lo, hi=hi, counts=cnt, mean=mu, var=var, p_low=cV / 20000)
    say('exponents per doubling from 2^%d to 2^%d: %s; mean %.2f (model e^gamma = %.2f), variance %.2f: variance/mean %.2f; a Poisson count is this even or more in %.1f%% of runs'
        % (lo, hi, cnt, mu, lam, var, var / mu, 100 * cV / 20000))

    # 2. remainders: same class mod m for consecutive exponents >= 61, against models A and B
    seq = [p for p in ex if p >= 61]
    pairs = len(seq) - 1
    real = {m: gaps_divisible(seq, m) for m in MODS}
    realw = {m: gaps_divisible(seq, m) for m in WHEELS}
    rng = random.Random(3)
    simA = {m: [] for m in MODS}
    simB = {m: [] for m in MODS}
    simBw = {m: [] for m in WHEELS}
    class4 = [0, 0]
    for s in range(sims):
        fa = [rand_prime(rng, p, False) for p in seq]
        fb = [rand_prime(rng, p, True) for p in seq]
        for m in MODS:
            simA[m].append(gaps_divisible(fa, m))
            simB[m].append(gaps_divisible(fb, m))
        for m in WHEELS:
            simBw[m].append(gaps_divisible(fb, m))
        class4[0] += sum(1 for q in fb if q % 4 == 1)
        class4[1] += len(fb)
    rows = []
    for m in MODS:
        a, b = simA[m], simB[m]
        ea, eb = sum(a) / sims, sum(b) / sims
        pb = min(1.0, 2 * min(sum(1 for x in b if x >= real[m]), sum(1 for x in b if x <= real[m])) / sims)
        rows.append(dict(m=m, real=real[m], A=ea, B=eb, pB=pb))
        say('  gap a multiple of %2d: %2d of %d; model A %.1f, model B %.1f (two-sided p = %.2f)' % (m, real[m], pairs, ea, eb, pb))
    out['mods'] = dict(pairs=pairs, start=61, rows=rows, sims=sims, K=K)
    say('exponents >= 61 that are 1 mod 4: %d of %d; model B gives %.1f%%' % (sum(1 for p in seq if p % 4 == 1), len(seq), 100 * class4[0] / class4[1]))
    out['class4'] = dict(real=sum(1 for p in seq if p % 4 == 1), n=len(seq), modelB=class4[0] / class4[1])
    # the gear version: every wheel m from 2 to 60
    wrows = []
    for m in WHEELS:
        b = simBw[m]
        eb = sum(b) / sims
        hi_ = sum(1 for x in b if x >= realw[m]) / sims
        lo_ = sum(1 for x in b if x <= realw[m]) / sims
        wrows.append(dict(m=m, real=realw[m], B=eb, p=min(1.0, 2 * min(hi_, lo_))))
    small = [w for w in wrows if w['p'] < 0.05]
    # how many wheels would fall below 0.05 in a sequence drawn from the model itself (leave-one-out)
    cnt_null = []
    for s in range(min(sims, 400)):
        c = 0
        for m in WHEELS:
            b = simBw[m]
            x = b[s]
            hi_ = (sum(1 for y in b if y >= x) - 1) / (sims - 1)
            lo_ = (sum(1 for y in b if y <= x) - 1) / (sims - 1)
            c += min(1.0, 2 * min(hi_, lo_)) < 0.05
        cnt_null.append(c)
    pglob = sum(1 for c in cnt_null if c >= len(small)) / len(cnt_null)
    out['wheels'] = dict(rows=wrows, below05=len(small), null_mean=sum(cnt_null) / len(cnt_null), p_global=pglob)
    say('wheels 2 to 60 (wheel m is home at gap g when m divides g): %d with p < 0.05 (%s); a model sequence gives %.1f on average; this many or more in %.0f%% of model sequences'
        % (len(small), ', '.join('%d: %d home, %.1f expected' % (w['m'], w['real'], w['B']) for w in small), sum(cnt_null) / len(cnt_null), 100 * pglob))

    # 3. the next exponent, if the model holds (heuristic)
    tested = 142351031
    q50 = tested * 2 ** (math.log(2) / lam)
    q90 = tested * 2 ** (math.log(10) / lam)
    out['next'] = dict(tested=tested, median=q50, q90=q90)
    say('if the model holds and no exponent below %d is missing: the next exponent has median %.3g and lies below %.3g with 90%%' % (tested, q50, q90))
    json.dump(out, open('mexp_gaps.json', 'w'))
    open('mexp_gaps.txt', 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
