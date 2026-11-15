"""Angle 2b: relations restricted to the tail k = 26..52 (a_k > 23209 for the real
sequence), where small-number coincidences are rare.
Ltail_D: number of k in 26..52 with a_k = c1 a_i + c2 a_j + d (i<j<k, 0<|c|<=8, |d|<=D)
         or a_k = c a_i + d (1<=c<=8).
Ltailtot_D: total number of such representations for k in 26..52."""
import json, sys, time
import numpy as np

CS = np.array([c for c in range(-8, 9) if c != 0], dtype=np.int64)
DS = (0, 1, 2, 10)


def tail_stats(a, kmin=25):
    a = np.array(a, dtype=np.int64)
    n = len(a)
    I, J = np.triu_indices(n, 1)
    keep = J < n - 1
    I, J = I[keep], J[keep]
    c1 = np.repeat(CS, CS.size)
    c2 = np.tile(CS, CS.size)
    V = a[I][:, None] * c1[None, :] + a[J][:, None] * c2[None, :]
    out = {}
    for D in DS:
        lo = np.searchsorted(a, V - D, side='left')
        hi = np.searchsorted(a, V + D, side='right')
        start = np.maximum(np.maximum(lo, (J + 1)[:, None]), kmin)
        cnt = np.clip(hi - start, 0, None)
        has = np.zeros(n, bool)
        for p_, c_ in zip(*np.nonzero(cnt)):
            has[start[p_, c_]:hi[p_, c_]] = True
        tot = int(cnt.sum())
        for c in range(1, 9):
            V1 = a[:-1] * c
            lo1 = np.searchsorted(a, V1 - D, side='left')
            hi1 = np.searchsorted(a, V1 + D, side='right')
            st1 = np.maximum(np.maximum(lo1, np.arange(1, n)), kmin)
            cn1 = np.clip(hi1 - st1, 0, None)
            tot += int(cn1.sum())
            for i in np.nonzero(cn1)[0]:
                has[st1[i]:hi1[i]] = True
        out['Ltail_D%d' % D] = int(has.sum())
        out['Ltailtot_D%d' % D] = tot
    return out


SIDES = {('Ltail_D%d' % D): 'hi' for D in DS}
SIDES.update({('Ltailtot_D%d' % D): 'hi' for D in DS})

if __name__ == '__main__':
    D = json.load(open('nulls.json'))
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    t0 = time.time()
    real = tail_stats(D['real'])
    print('real', real)
    out = {s: {'side': SIDES[s], 'real': real[s], 'a': [], 'b': []} for s in SIDES}
    for fam in ('a', 'b'):
        for seq in D[fam][:nmax]:
            st = tail_stats(seq)
            for s in SIDES:
                out[s][fam].append(st[s])
    json.dump(out, open('angle2b.json', 'w'))
    print('done', time.time() - t0)
