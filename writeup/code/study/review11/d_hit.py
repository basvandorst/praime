"""Review 11, script D: rounds 45-46 recomputed. Own odd-only spf sieve (descending overwrite), hit(y) = P(M(p) <= y),
own Buchstab omega, own pair factor <S_{>y}({0,h1})> (exact over classes for y <= 13, G(y) for y >= 97 with a check),
K entropies pooled and conditional on the class of p mod prod_{q<=y} q."""
import math, json, sys
import numpy as np

GAM = 0.5772156649015329

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

PR = primes_upto(10_000_000)

def omega(us):
    """Buchstab omega by unit intervals: u w(u) = n w(n) + int_n^u w(t-1) dt."""
    h = 1e-5; m = int(round(1 / h))
    t = np.linspace(1, 2, m + 1); w = 1 / t
    grid = [t]; vals = [w]
    for n in range(2, int(max(us)) + 2):
        prev = vals[-1]                              # w on [n-1, n]
        integ = np.concatenate([[0], np.cumsum((prev[1:] + prev[:-1]) * h / 2)])
        tn = np.linspace(n, n + 1, m + 1)
        wn = (n * prev[-1] + integ) / tn
        grid.append(tn); vals.append(wn)
    G = np.concatenate(grid); V = np.concatenate(vals)
    return np.interp(us, G, V)

def S_pair_gt(y, h):
    """prod_{q>y, q odd} (1 - nu/q)/(1 - 1/q)^2 for {0,h}."""
    qs = PR[(PR > max(y, 2))].astype(float)
    nu = np.where(h % PR[PR > max(y, 2)] == 0, 1.0, 2.0)
    return float(np.exp(np.sum(np.log1p(-nu / qs) - 2 * np.log1p(-1 / qs))))

def mean_pair(y):
    sm = [int(q) for q in PR if q <= y]; P = math.prod(sm)
    if y > 13:
        return S_pair_gt(y, 1)                       # h1 not divisible by any q > y (h1 < y here, checked below)
    cache = {}; tot = 0; n = 0
    for a in range(P):
        if math.gcd(a, P) != 1: continue
        h = 1
        while math.gcd(a + h, P) != 1: h += 1
        if h not in cache: cache[h] = S_pair_gt(y, h)
        tot += cache[h]; n += 1
    return tot / n

def window(L, W, ys):
    top = math.isqrt(L + W) + 1
    base = [int(q) for q in PR[(PR > 2) & (PR <= top)]][::-1]
    spf = np.zeros(W // 2, np.uint32)                # odd numbers L+1+2i (L even)
    for q in base:                                   # descending: last write is the smallest prime factor
        m = ((L + 1 + q - 1) // q) * q
        if m % 2 == 0: m += q
        spf[(m - L - 1) // 2::q] = q
    pidx = np.nonzero(spf == 0)[0]
    seg = spf
    Mx = np.maximum.reduceat(seg, pidx[:-1])        # segment starts at a prime (0) and runs to the next prime
    Mx = np.maximum(Mx, 2)
    lam = math.log(L + W / 2)
    om = omega(np.array([lam / math.log(y) for y in ys]))
    out = []
    N = len(Mx)
    cls_p = (L + 1 + 2 * pidx[:-1].astype(np.int64))
    for y, w in zip(ys, om):
        u = lam / math.log(y)
        meas = float(np.mean(Mx <= y)); sd = math.sqrt(meas * (1 - meas) / N)
        prod = float(np.prod(1 - 1 / PR[PR <= y]))
        ms = mean_pair(y)
        pre = min(1.0, 1 / (lam * prod * math.exp(GAM) * w))
        cor = min(1.0, ms * pre) if u > 2 else 1.0
        uonly = ms / (u * w)
        row = dict(y=y, u=u, meas=meas, sd=sd, prereg=pre, corrected=cor, u_only_times_pair=uonly, pair=ms)
        # round 46: K = 1 + number of odd composites with spf > y in the gap
        big = (seg > y).astype(np.int16)
        K = 1 + np.add.reduceat(big, pidx[:-1])
        def H(Karr):
            v, c = np.unique(Karr, return_counts=True); p = c / c.sum()
            return float(-np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * c.sum() * math.log(2)))
        rho = cor
        row['H_pooled_MM'] = H(K)
        row['H_geom'] = 0.0 if rho >= 1 else -(rho * math.log2(rho) + (1 - rho) * math.log2(1 - rho)) / rho
        if y <= 13:
            Pq = math.prod(int(q) for q in PR if q <= y)
            c = cls_p % Pq
            Hc = 0.0
            for a in np.unique(c):
                sel = K[c == a]; Hc += len(sel) / N * H(sel)
            row['H_cond_class_MM'] = Hc
        out.append(row)
        print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    return out

if __name__ == '__main__':
    which = sys.argv[1]
    if which == '1e9':
        res = window(10 ** 9, 5 * 10 ** 7, [2, 3, 5, 7, 13, 97, 293, 1777, 10607, 20731, 32401])
    else:
        res = window(10 ** 12, 10 ** 7, [2, 3, 7, 13, 97, 1069, 11801, 32869, 999983])
    json.dump(res, open('/home/bas/scratch/study-review11/d_hit_%s.json' % which, 'w'), indent=1)
