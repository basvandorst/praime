"""Part 2: P-1 on top of trial factoring.  Reads pm1_<range>.json (pm1.py) and the kmin files.
Validation of every P-1 factor g (plain integers + sympy): 1 < g < M_p and 2^p = 1 mod g (so M_p is proved composite);
if g is prime: g = 2kp + 1 and k factors as the bounds say (stage 1 at B1: k is B1-powersmooth up to the factor 2;
stage 2 at B2: k / (one prime r <= B2) is so).
Cost in LL equivalents = operations / (p - 2): stage 1 = bitlen(E_B1) squarings, stage 2 = 2 multiplications per prime
in (B1, B2] (a multiplication counted as one squaring), stage 2 skipped when stage 1 found a factor."""
import json, sys, math, numpy as np, gmpy2
from sympy import factorint, isprime
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243]
E = set(MEXP)

def smooth_ok(k, B1, B2=None):
    f = factorint(k)
    big = [(l, e) for l, e in f.items() if l ** e > B1 and not (l == 2)]
    if B2 is None: return not big
    return len(big) == 0 or (len(big) == 1 and big[0][1] == 1 and big[0][0] <= B2)

def report(jf, kf, Ks, X, pmin=0):
    d = json.load(open(jf)); res = {int(p): r for p, r in d['res'].items()}
    S1, S2 = {int(b): v for b, v in d['S1COST'].items()}, {int(c): v for c, v in d['S2COST'].items()}
    lev, S2B1, CP, Kd = d['levels'], d['S2B1'], d['checkpoints'], d['K']
    P, km = np.load(kf); P = [int(p) for p in P]; km = [int(k) for k in km]
    kmin = dict(zip(P, km))
    nsm = nbad = 0
    for p, r in res.items():                                   # validation
        if r['g'] is None: continue
        g = int(r['g']); M = (1 << p) - 1
        assert 1 < g < M and pow(2, p, g) == 1, p
        assert (g - 1) % (2 * p) == 0 or not isprime(g), p
        if isprime(g):
            k = (g - 1) // (2 * p)
            assert k > Kd, (p, k)
            ok = smooth_ok(k, r['s1']) if r['s1'] is not None else smooth_ok(k, S2B1, r['s2'])
            nsm += 1; nbad += (not ok)
    print(f'{jf}: {sum(1 for r in res.values() if r["g"])} factors found, all divide M_p; {nsm} prime factors, k smooth as the bounds say for all but {nbad} (those need ord_q(3) < q - 1)')
    print(f'  p in [{max(pmin, d["lo"])},{d["hi"]}), TF depth K, P-1 config: survivors -> LL tests, extra removed by P-1, hits, LL/exp, P-1 cost (LL-eq), cost per removal (LL-eq)')
    for K in Ks:
        base = [p for p in res if (kmin[p] == 0 or kmin[p] > K) and p >= pmin]
        hits = sum(p in E for p in base)
        assert K >= Kd and (d.get('stride', 1) > 1 or all(p in res for p in P if d['lo'] <= p < d['hi'] and (kmin[p] == 0 or kmin[p] > K)))
        if d.get('stride', 1) > 1: print(f'  SAMPLE: every {d["stride"]}-th survivor only')
        configs = [('none', None, None)] + [(f'B1={B}', B, None) for B in lev] + [(f'B1={S2B1},B2={c}', S2B1, c) for c in CP]
        for name, B1, B2 in configs:
            if B1 is None:
                rem = 0; cost = 0.0
            else:
                rem = 0; cost = 0.0
                for p in base:
                    r = res[p]
                    f1 = r['s1'] is not None and r['s1'] <= B1
                    f2 = B2 is not None and not f1 and r['s2'] is not None and r['s2'] <= B2
                    c = S1[B1] + (S2[B2] if (B2 is not None and not f1) else 0)
                    cost += c / (p - 2) if p > 2 else 0
                    rem += f1 or f2
            n = len(base) - rem + (1 if max(pmin, d['lo']) <= 3 else 0)        # + p = 2 for the range starting at 2/3
            h = hits + (1 if max(pmin, d['lo']) <= 3 else 0)
            per = f'{cost / rem:8.2f}' if rem else '       -'
            print(f'  K={K:>8} {name:<22} LL {n:>5}  removed {rem:>4}  hits {h}  LL/exp {(n / h if h else float('inf')):7.1f}  P-1 {cost:9.1f}  per removal {per}')

if __name__ == '__main__':
    if sys.argv[1] == '1e4':
        report('pm1_1e4.json', 'kmin_10000_10000000.npy', [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7], 10000)
        report('pm1_1e4.json', 'kmin_10000_10000000.npy', [10 ** 4, 10 ** 6], 10000, pmin=1000)
    try:
        if sys.argv[1] == '1e5':
            for c in 'abc': report(f'pm1_1e5{c}.json', 'kmin_100000_1000000.npy', [10 ** 6], 100000)
    except FileNotFoundError: pass
