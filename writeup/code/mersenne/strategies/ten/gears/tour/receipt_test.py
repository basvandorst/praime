import json, random, sys, time
import gmpy2
from gmpy2 import mpz
from sympy import primerange, isprime
from receipt import One, Two, run, build, check, forge, power, lucas_lehmer, halves

rng = random.Random(1811)
out = {}
lines = []


def say(s):
    print(s, flush=True)
    lines.append(s)


TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
t0 = time.time()
agree = {'one': 0, 'two': 0}
acc = {'one': 0, 'two': 0}
tot = 0
exps = []
for p in primerange(3, TOP + 1):
    ll = lucas_lehmer(p)
    tot += 1
    if ll:
        exps.append(p)
    for cls, key in ((One, 'one'), (Two, 'two')):
        G = cls(p)
        k = 3
        y, saved, top = run(G, k)
        if G.verdict(y) == ll:
            agree[key] += 1
        if check(G, y, build(G, k, y, saved)):
            acc[key] += 1
say('odd primes p up to %d: %d; verdict equal to Lucas-Lehmer: one cell %d, two cells %d; honest receipts accepted: one cell %d, two cells %d; exponents found: %s (%.0f s)' % (TOP, tot, agree['one'], agree['two'], acc['one'], acc['two'], exps, time.time() - t0))
out['agree'] = dict(top=TOP, primes=tot, one=agree['one'], two=agree['two'], acc_one=acc['one'], acc_two=acc['two'], exps=exps)

comp_pass = [n for n in range(4, 2000) if not isprime(n) and One(n).verdict(run(One(n), 1)[0])]
say('composite exponents n from 4 to 1999 whose 2^n - 1 passes the one-cell test: %s' % comp_pass)
out['composite_n_pass'] = comp_pass

rows = []
for p in [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else '9689,21701,44497,86243,110503').split(',')]:
    for cls, name in ((One, 'one cell'), (Two, 'two cells')):
        for k in (6, 8):
            G = cls(p)
            t = time.time(); y, saved, top = run(G, k); t1 = time.time() - t
            cnt = [0]
            t = time.time(); mus = build(G, k, y, saved, count=cnt); t2 = time.time() - t
            c2 = [0]
            t = time.time(); ok = check(G, y, mus, count=c2); t3 = time.time() - t
            assert ok and G.verdict(y)
            row = dict(p=p, group=name, k=k, squarings=top, saved=len(saved), receipt_products=cnt[0], check_squarings=c2[0], test_s=round(t1, 3), receipt_s=round(t2, 3), check_s=round(t3, 4), ratio=round(t1 / t3, 1))
            rows.append(row)
            say('p = %d, %s, %d halvings: test %d squarings %.2f s; receipt: %d numbers, built from %d saved residues in %.2f s; check: %d squarings and %d small powers, %.3f s; the check is %.0f times faster than the test' % (p, name, k, top, t1, k, len(saved), t2, c2[0], 2 * k, t3, t1 / t3))
out['timing'] = rows

p = 9689
res = {}
for cls, key in ((One, 'one'), (Two, 'two')):
    G = cls(p)
    k = 6
    y, saved, top = run(G, k)
    mus = build(G, k, y, saved)
    assert check(G, y, mus)
    bump = (lambda v: (v + 1) % G.N) if key == 'one' else (lambda v: ((v[0] + 1) % G.N, v[1]))
    rnd = (lambda: mpz(rng.randrange(2, int(G.N)))) if key == 'one' else (lambda: (mpz(rng.randrange(2, int(G.N))), mpz(rng.randrange(2, int(G.N)))))
    wrong_y = sum(1 for _ in range(50) if check(G, rnd(), mus)) + (1 if check(G, bump(y), mus) else 0)
    wrong_mu = 0
    for i in range(k):
        m2 = list(mus); m2[i] = bump(m2[i])
        wrong_mu += 1 if check(G, y, m2) else 0
        m2 = list(mus); m2[i] = rnd()
        wrong_mu += 1 if check(G, y, m2) else 0
    flips = 0
    same = 0
    NF = 100
    for _ in range(NF):
        st = rng.randrange(1, G.T + 1)
        y2, saved2, _ = run(G, k, flip=(st, rng.randrange(p)))
        if y2 == y:
            same += 1
        if check(G, y2, build(G, k, y2, saved2)):
            flips += 1
    res[key] = dict(wrong_y_accepted=wrong_y, wrong_y_tried=51, wrong_mu_accepted=wrong_mu, wrong_mu_tried=2 * k, flips=NF, flips_accepted=flips, flips_same_result=same)
    say('p = %d, %s: wrong result accepted %d of 51; one changed receipt number accepted %d of %d; a flipped bit during the test (%d runs): receipt accepted %d times, result unchanged %d times' % (p, key, wrong_y, wrong_mu, 2 * k, NF, flips, same))
out['tamper'] = res

forg = []
for p in (521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689):
    G = One(p)
    k = 4
    y, saved, top = run(G, k)
    N = G.N
    g = 2
    while True:
        g += 1
        w = gmpy2.powmod(g, (N - 1) // 3, N)
        if w != 1:
            break
    for name, eps, d in (('-1 (order 2)', N - 1, 2), ('a cube root of 1 (order 3)', w, 3)):
        tried = hit = accepted = flips_verdict = 0
        for nonce in range(40):
            tried += 1
            fz = forge(G, k, y, saved, eps, d, nonce)
            if fz is None:
                continue
            hit += 1
            ys, mus = fz
            assert ys != y
            if check(G, ys, mus, nonce):
                accepted += 1
                if G.verdict(ys) != G.verdict(y):
                    flips_verdict += 1
        forg.append(dict(p=p, group='one', factor=name, tried=tried, built=hit, accepted=accepted, verdict_changed=flips_verdict))
    G = Two(p)
    y, saved, top = run(G, k)
    for name, eps, d in (('-1 (order 2)', (N - 1, mpz(0)), 2), ('the true result itself (order 4)', y, 4)):
        tried = hit = accepted = flips_verdict = 0
        for nonce in range(40):
            tried += 1
            fz = forge(G, k, y, saved, eps, d, nonce)
            if fz is None:
                continue
            hit += 1
            ys, mus = fz
            assert ys != y
            if check(G, ys, mus, nonce):
                accepted += 1
                if G.verdict(ys) != G.verdict(y):
                    flips_verdict += 1
        forg.append(dict(p=p, group='two', factor=name, tried=tried, built=hit, accepted=accepted, verdict_changed=flips_verdict))
for grp in ('one', 'two'):
    for name in sorted(set(r['factor'] for r in forg if r['group'] == grp)):
        rr = [r for r in forg if r['group'] == grp and r['factor'] == name]
        say('forgery on prime 2^p - 1, %s cell(s), true result times %s: %d attempts on %d exponents; a false receipt was built in %d and the checker accepted %d of them; the verdict changed from prime to not prime in %d' % (
            grp, name, sum(r['tried'] for r in rr), len(rr), sum(r['built'] for r in rr), sum(r['accepted'] for r in rr), sum(r['verdict_changed'] for r in rr)))
out['forgery'] = forg

S = 20
safe = dict(one_forged_flagged=0, one_forged=0, two_forged_flagged=0, two_forged=0)
for p in (521, 1279, 4423):
    N = (mpz(1) << p) - 1
    L = 2
    m = (N - 1) // 2
    for q in primerange(3, 1 << 16):
        while m % q == 0:
            m //= q
            L *= q
    G = One(p)
    y, saved, _ = run(G, 4)
    g = 2
    while True:
        g += 1
        w = gmpy2.powmod(g, (N - 1) // 3, N)
        if w != 1:
            break
    for eps in (N - 1, w, w * w % N, (N - 1) * w % N):
        ys = y * eps % N
        safe['one_forged'] += 1
        if gmpy2.powmod(ys, L, N) == gmpy2.powmod(9, L, N):
            safe['one_forged_flagged'] += 1
    G = Two(p)
    y, saved, _ = run(G, 4)
    for eps in ((N - 1, mpz(0)), y, G.mul(y, y), G.mul(G.mul(y, y), y)):
        ys = G.mul(y, eps)
        v = ys
        for _ in range(S):
            v = G.sq(v)
        safe['two_forged'] += 1
        if v == G.one():
            safe['two_forged_flagged'] += 1
flag_comp = {'one': 0, 'two': 0}
ncomp = 0
for p in primerange(3, 2000):
    if lucas_lehmer(p):
        continue
    ncomp += 1
    N = (mpz(1) << p) - 1
    L = 2
    m = (N - 1) // 2
    for q in primerange(3, 1 << 12):
        while m % q == 0:
            m //= q
            L *= q
    y = run(One(p), 1)[0]
    if gmpy2.powmod(y, L, N) == gmpy2.powmod(9, L, N):
        flag_comp['one'] += 1
    G = Two(p)
    v = run(G, 1)[0]
    for _ in range(S):
        v = G.sq(v)
    if v == G.one():
        flag_comp['two'] += 1
say('safer verdict (compare after raising to a power that kills factors of small order): forged results flagged for a full test: one cell %d of %d, two cells %d of %d; composite 2^p - 1 with prime p < 2000 flagged without need: one cell %d of %d, two cells %d of %d' % (
    safe['one_forged_flagged'], safe['one_forged'], safe['two_forged_flagged'], safe['two_forged'], flag_comp['one'], ncomp, flag_comp['two'], ncomp))
out['safe'] = dict(safe, composites=ncomp, flagged_one=flag_comp['one'], flagged_two=flag_comp['two'], S=S)
json.dump(out, open('receipt_test.json', 'w'))
open('receipt_test.txt', 'w').write('\n'.join(lines) + '\n')
