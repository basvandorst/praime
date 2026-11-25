"""genus-2 candidates for the degree-4 follow-up: box y^2 + h y = f, f in {-1,0,1}^7 (deg 5 or 6), h one of the 16 0/1-polynomials
of degree <= 3 (as round 151); genus2red conductor N_red < 4000 (a lower bound when the 2-part is unknown); one representative per
fingerprint (N_red, #C(F_p) data via hyperellcharpoly at the six smallest good odd primes); then lfungenus2: conductor N_L from
lfunparams, lfuncheckfeq <= -20, N_L <= 3000, one per (N_L, a_1..a_60), primitive iff mean |lambda_p|^2 over good p <= 3000
in [0.75, 1.25]; excluded: every genus-2 L-function of round 1 and round 151 (tests, proxies, 169, 249, 277, 529, 961).
Writes genus2.json (all valid rows) for the seeded draw."""
import os, json, time, itertools, math
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
t0 = time.time()
hs = ['+'.join(['x^%d' % j for j in range(4) if bits >> j & 1]) or '0' for bits in range(16)]
cand = {}
n = 0
for coefs in itertools.product((-1, 0, 1), repeat=7):
    if coefs[6] == 0 and coefs[5] == 0: continue
    f = '+'.join('(%d)*x^%d' % (c, j) for j, c in enumerate(coefs) if c != 0)
    for h in hs:
        try:
            F = pari('4*(%s)+(%s)^2' % (f, h))
            if int(pari.poldegree(F)) < 5 or int(pari.poldegree(pari.gcd(F, pari.deriv(F)))) > 0: continue
            N = int(pari('genus2red([%s,%s])' % (f, h))[0])
        except Exception:
            continue
        n += 1
        if N >= 4000: continue
        good = [p for p in (3, 5, 7, 11, 13, 17, 19, 23, 29) if N % p][:6]
        fp = (N,) + tuple(str(pari('hyperellcharpoly(Mod(1,%d)*(4*(%s)+(%s)^2))' % (p, f, h))) for p in good)
        if fp not in cand: cand[fp] = (N, f, h)
print('curves', n, 'fingerprints with N_red < 4000:', len(cand), 'secs', round(time.time() - t0, 1), flush=True)
R151 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'blind', 'testset.json')))
old = [[r[1], r[2]] for r in R151['T4b']['test'] + R151['T4b']['proxy']] + [['x^2+x', 'x^3+1'], ['x^5-2*x^3+2*x^2-x', '1'], ['x^5+x^4', 'x^3+x+1']]
old += [[r[1], r[2]] for r in R151['T4b']['all_valid'] if r[0] in (529, 961)]
oldkeys = set()
for f, h in old:
    L = pari('lfungenus2([%s,%s])' % (f, h)); an = pari.lfunan(L, 60)
    oldkeys.add((int(pari.lfunparams(L)[0]),) + tuple(int(x) for x in an))
rows = []; seen = set()
for fp, (N, f, h) in sorted(cand.items()):
    try:
        L = pari('lfungenus2([%s,%s])' % (f, h)); NL = int(pari.lfunparams(L)[0])
        if NL > 3000: continue
        an = pari.lfunan(L, 3000)
        key = (NL,) + tuple(int(an[i]) for i in range(60))
        if key in seen or key in oldkeys: continue
        feq = float(pari.lfuncheckfeq(L))
    except Exception:
        continue
    seen.add(key)
    if feq > -20: continue
    ps = [int(p) for p in pari.primes([2, 3000]) if NL % int(p)]
    m2 = float(np.mean([(float(an[p - 1]) / math.sqrt(p)) ** 2 for p in ps]))
    rows.append([NL, f, h, round(m2, 4), round(feq)])
rows.sort()
prim = [r for r in rows if 0.75 <= r[3] <= 1.25]
json.dump(dict(valid=rows, primitive=prim), open('genus2.json', 'w'), indent=1)
print('valid L-functions with N_L <= 3000:', len(rows), ' primitive:', len(prim), ' secs', round(time.time() - t0, 1))
print('conductors (primitive):', [r[0] for r in prim])
