r"""Angle 2: relations among terms.

L2(D): a_k = c1*a_i + c2*a_j + d, i < j < k, c1, c2 in [-8, 8] \ {0}, |d| <= D.
L1(D): a_k = c*a_i + d, i < k, c in [1, 8], |d| <= D.
  stats per D in {0,1,2,10}: total number of representations, number of k >= 9
  (a_k beyond the 8 smallest terms) with at least one representation.
Out of sample (lag patterns): a_k = c1*a_{k-r} + c2*a_{k-s} + d, 1 <= r < s <= 6,
  c1, c2 in [-8,8]\{0}, |d| <= D; the pattern with most hits on k <= 25 is chosen,
  ties broken by order of enumeration; stat = hits for k = 26..52.
MUL: a_k = round(a_i*a_j / a_l), i <= j, l not in {i, j}, all indices < k.
POW: a_k = round(a_i^(p/q)), q = 1..4, p/q > 1, i < k (exact), and the looser
  |a_i^(p/q) - a_k| / a_k < 1e-3.
"""
import json, sys, time
import numpy as np

CS = np.array([c for c in range(-8, 9) if c != 0], dtype=np.int64)
DS = (0, 1, 2, 10)


def lin_stats(a, want_list=False):
    a = np.array(a, dtype=np.int64)
    n = len(a)
    I, J = np.triu_indices(n, 1)  # i < j
    keep = J < n - 1
    I, J = I[keep], J[keep]
    out = {}
    lists = {}
    c1 = np.repeat(CS, CS.size)
    c2 = np.tile(CS, CS.size)
    V = (a[I][:, None] * c1[None, :] + a[J][:, None] * c2[None, :])  # (pairs, 256)
    for D in DS:
        lo = np.searchsorted(a, V - D, side='left')
        hi = np.searchsorted(a, V + D, side='right')
        start = np.maximum(lo, (J + 1)[:, None])
        cnt = np.clip(hi - start, 0, None)
        tot = int(cnt.sum())
        hasrep = np.zeros(n, dtype=bool)
        recs = []
        pi, ci = np.nonzero(cnt)
        for p_, c_ in zip(pi, ci):
            for kk in range(start[p_, c_], hi[p_, c_]):
                hasrep[kk] = True
                if want_list:
                    recs.append((int(kk), int(I[p_]), int(J[p_]), int(c1[c_]), int(c2[c_]),
                                 int(a[kk] - V[p_, c_])))
        # single-term
        tot1 = 0
        for c in range(1, 9):
            V1 = a[:-1] * c
            lo1 = np.searchsorted(a, V1 - D, side='left')
            hi1 = np.searchsorted(a, V1 + D, side='right')
            st1 = np.maximum(lo1, np.arange(1, n))
            cn1 = np.clip(hi1 - st1, 0, None)
            tot1 += int(cn1.sum())
            for i in np.nonzero(cn1)[0]:
                for kk in range(st1[i], hi1[i]):
                    hasrep[kk] = True
                    if want_list:
                        recs.append((int(kk), int(i), -1, c, 0, int(a[kk] - V1[i])))
        out['L2tot_D%d' % D] = tot
        out['L1tot_D%d' % D] = tot1
        out['Lhas_D%d' % D] = int(hasrep[8:].sum())
        if want_list:
            lists[D] = recs
    return out, lists


def lag_oos(a, D, n_in=25, smax=6):
    a = np.array(a, dtype=np.int64)
    n = len(a)
    best = None
    for r in range(1, smax):
        for s in range(r + 1, smax + 1):
            k = np.arange(s, n)  # 0-based k with k-s >= 0
            V = a[k - r][:, None] * np.repeat(CS, CS.size)[None, :] + a[k - s][:, None] * np.tile(CS, CS.size)[None, :]
            hit = np.abs(a[k][:, None] - V) <= D
            ins = hit[k < n_in].sum(0)
            oos = hit[k >= n_in].sum(0)
            j = int(np.argmax(ins))
            if best is None or ins[j] > best[0]:
                best = (int(ins[j]), int(oos[j]), r, s, int(np.repeat(CS, CS.size)[j]), int(np.tile(CS, CS.size)[j]))
    return best


def mul_stats(a):
    a = np.array(a, dtype=np.int64)
    n = len(a)
    I, J = np.triu_indices(n - 1, 0)  # i <= j < n-1
    tot = 0
    hasrep = np.zeros(n, dtype=bool)
    prods = a[I] * a[J]  # < 2^63 since a < 2^28
    for l in range(n - 1):
        ok = (I != l) & (J != l)
        x = prods[ok]
        y = int(a[l])
        v = (2 * x + y) // (2 * y)
        m = np.maximum(np.maximum(I[ok], J[ok]), l)
        pos = np.searchsorted(a, v)
        pos_c = np.minimum(pos, n - 1)
        hit = (a[pos_c] == v) & (pos_c > m)
        tot += int(hit.sum())
        hasrep[pos_c[hit]] = True
    return tot, int(hasrep[8:].sum())


FR = sorted({(p, q) for q in range(1, 5) for p in range(q + 1, 8 * q + 1) if np.gcd(p, q) == 1})


def pow_stats(a):
    a = np.array(a, dtype=np.float64)
    ai = np.array(a, dtype=np.int64)
    n = len(a)
    exact = 0
    loose = 0
    upper = np.triu(np.ones((n, n), dtype=bool), 1)  # [i, k] with i < k
    for p, q in FR:
        v = a ** (p / q)
        vr = np.round(v)
        exact += int(((vr[:, None] == a[None, :]) & upper).sum())
        loose += int(((np.abs(v[:, None] - a[None, :]) / a[None, :] < 1e-3) & upper).sum())
    return exact, loose


def all_stats(a):
    out, _ = lin_stats(a)
    for D in DS:
        b = lag_oos(a, D)
        out['lag_ins_D%d' % D] = b[0]
        out['lag_oos_D%d' % D] = b[1]
    out['mul_tot'], out['mul_has'] = mul_stats(a)
    out['pow_exact'], out['pow_loose'] = pow_stats(a)
    return out


def sides():
    s = {}
    for D in DS:
        for k in ('L2tot_D%d', 'L1tot_D%d', 'Lhas_D%d', 'lag_ins_D%d', 'lag_oos_D%d'):
            s[k % D] = 'hi'
    for k in ('mul_tot', 'mul_has', 'pow_exact', 'pow_loose'):
        s[k] = 'hi'
    return s


if __name__ == '__main__':
    D = json.load(open('nulls.json'))
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    t0 = time.time()
    A = D['real']
    real = all_stats(A)
    SIDES = sides()
    with open('angle2_real.txt', 'w') as rep:
        rep.write('real stats: %s\n' % real)
        _, lists = lin_stats(A, want_list=True)
        for Dv in (0, 1, 2):
            rep.write('\nD=%d: representations a_k = c1 a_i + c2 a_j + d (k >= 13, i.e. a_k >= 521)\n' % Dv)
            for (kk, i, j, c1, c2, d) in sorted(lists[Dv]):
                if kk >= 12:
                    if j >= 0:
                        rep.write('  %d = %d*%d + %d*%d + %d\n' % (A[kk], c1, A[i], c2, A[j], d))
                    else:
                        rep.write('  %d = %d*%d + %d\n' % (A[kk], c1, A[i], d))
            per_k = {}
            for (kk, *_r) in lists[Dv]:
                per_k[A[kk]] = per_k.get(A[kk], 0) + 1
            rep.write('  count per a_k: %s\n' % per_k)
        for Dv in DS:
            rep.write('lag pattern D=%d: (ins hits, oos hits, r, s, c1, c2) = %s\n' % (Dv, lag_oos(A, Dv)))
    out = {s: {'side': SIDES[s], 'real': real[s], 'a': [], 'b': []} for s in SIDES}
    for fam in ('a', 'b'):
        for seq in D[fam][:nmax]:
            st = all_stats(seq)
            for s in SIDES:
                out[s][fam].append(st[s])
    json.dump(out, open('angle2.json', 'w'))
    print('done', time.time() - t0)
