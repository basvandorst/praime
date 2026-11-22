"""Independent data side of rounds 19, 21, 22: own segmented sieve, sliding counts, gaps, block t9 against own exact law."""
import math, json, sys
import numpy as np

OUT = '/home/bas/scratch/study-review3/'
E = json.load(open(OUT + 'enum.json'))


def base_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


def sieve(L, W, seg=10_000_000):
    """bool mask of primes in [L, L+W)."""
    bp = base_primes(int(math.isqrt(L + W)) + 1)
    M = np.ones(W, dtype=bool)
    for lo in range(L, L + W, seg):
        hi = min(lo + seg, L + W)
        s = M[lo - L:hi - L]
        for p in bp:
            p = int(p)
            if p * p >= hi:
                break
            st = max(p * p, (lo + p - 1) // p * p)
            s[st - lo::p] = False
    if L < 2:
        M[:2 - L] = False
    return M


def pk(T, lam, K):
    return np.array([sum(math.comb(j, k) * (-1) ** (j - k) * T[j] / lam ** j for j in range(k, len(T))) for k in range(K + 1)])


def share(U, lam):
    return sum((-1) ** j * U[j] / lam ** (j + 1) for j in range(len(U)))


def tstat(r):
    r = np.asarray(r)
    return float(r.mean() / (r.std(ddof=1) / math.sqrt(len(r))))


def simpson(f, a, b, n=200):
    x = np.linspace(a, b, n + 1); y = np.array([f(t) for t in x])
    return float((b - a) / n / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum()))


def window(name, L, W, do_gaps):
    M = sieve(L, W + 100)
    res = {}
    for H, K in ((30, 7), (60, 10)):
        T = E['H%d' % H]['T']
        emp, th = [], []
        for b in range(10):
            lo, n = b * W // 10, W // 10
            c = np.concatenate([[0], np.cumsum(M[lo:lo + n + H - 1], dtype=np.int32)])
            cnt = c[H:] - c[:-H]                             # primes in [L+lo+i, L+lo+i+H), i < n
            assert len(cnt) == n
            emp.append(np.bincount(cnt, minlength=K + 1)[:K + 1] / n)
            a, bb = L + lo, L + lo + n
            th.append(np.array([simpson(lambda t: pk(T, math.log(t), K)[k], a, bb, 40) for k in range(K + 1)]) / (bb - a))
        emp = np.array(emp); th = np.array(th)
        t9 = [tstat(emp[:, k] - th[:, k]) if th[:, k].mean() > 1e-6 else None for k in range(K + 1)]
        res['H%d' % H] = {'data': emp.mean(0).tolist(), 'law': th.mean(0).tolist(), 't9': t9}
        print(name, H, ' '.join('%d:%.5f/%.5f(%s)' % (k, emp.mean(0)[k], th.mean(0)[k], '%.2f' % t9[k] if t9[k] is not None else '-')
                                for k in range(K + 1)), flush=True)
    if do_gaps:
        P = np.nonzero(M[:W])[0].astype(np.int64) + L
        P = np.append(P, L + W + np.nonzero(M[W:])[0][0])   # include the closing gap (differs from r22 by one gap)
        g = np.diff(P)
        U = {int(h): v for h, v in E['U'].items()}
        st = {}
        for b in range(10):
            a, bb = L + b * (W // 10), L + (b + 1) * (W // 10)
            sel = (P[:-1] >= a) & (P[:-1] < bb); gb = g[sel]
            dens = simpson(lambda t: 1 / math.log(t), a, bb, 40)
            sh = {h: simpson(lambda t: share(U[h], math.log(t)) / math.log(t), a, bb, 40) / dens for h in U}
            for h in U:
                st.setdefault('gap%d' % h, []).append((np.mean(gb == h), sh[h]))
            st.setdefault('tail54', []).append((np.mean(gb >= 54), 1 - sum(v for h, v in sh.items() if h <= 52)))
            st.setdefault('tail72', []).append((np.mean(gb >= 72), 1 - sum(sh.values())))
        out = {}
        for k, v in st.items():
            v = np.array(v); out[k] = {'data': float(v[:, 0].mean()), 'law': float(v[:, 1].mean()), 't9': tstat(v[:, 0] - v[:, 1])}
        res['gaps'] = out
        print(name, 'tail54', out['tail54'], 'tail72', out['tail72'], flush=True)
        print(name, 'gap misses', [(k, round(d['t9'], 2)) for k, d in out.items() if abs(d['t9']) > 3],
              'max|t9|', max(abs(d['t9']) for d in out.values()), flush=True)
    return res


if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'train':
        r = window('train', 10 ** 7, 9 * 10 ** 7, True)
    elif which == 'A':
        r = window('A', 10 ** 10, 5 * 10 ** 7, True)
    elif which == 'jc':
        M = sieve(0, 10 ** 8 + 1000)
        P = np.nonzero(M)[0]
        g = np.diff(P); last = P[1:] <= 10 ** 8              # pairs with the larger prime <= 1e8
        n6 = int(np.sum((g == 6) & last)); n30 = int(np.sum((g == 30) & last))
        n6b = int(np.sum((g == 6) & (P[:-1] <= 10 ** 8))); n30b = int(np.sum((g == 30) & (P[:-1] <= 10 ** 8)))
        print('N(6,1e8)', n6, 'N(30,1e8)', n30, 'ratio', n30 / n6, '(by smaller prime:', n6b, n30b, ')')
        r = {'N6': n6, 'N30': n30}
    json.dump(r, open(OUT + 'data_%s.json' % which, 'w'), indent=1)
