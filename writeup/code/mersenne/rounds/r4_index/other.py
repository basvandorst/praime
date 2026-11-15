"""Closed forms, linear recurrences and relations on every derived sequence (real row 0 and the 300 null rows).

Closed forms b_n = g(n) (strip arrangement of formulas.py, section 6), n = 1..len:
  F1 floor(c r^n), F3 floor(c n^k), F4 round(c r^n): best exact hits (any two terms can always be hit);
  OOS: fits chosen on the in-sample terms (max in-sample hits), max hits among the other terms over the tied fits.
Linear recurrences (relations/angle1_linrec.py):
  bm_zero : #m >= 8 where the shortest recurrence of b_1..b_m (Berlekamp-Massey mod 2^61-1) predicts b_{m+1};
  bm_def  : max_m (ceil(m/2) - L(m));
  si_*    : b_n = sum_{i<=k} c_i b_{n-i} (+ d), k <= 3, |c_i| <= 10: longest run, run from the start,
            OOS hits of the recurrence with the longest run ending at the last in-sample term.
Relations (relations/angle2_relations.py, generalised to unsorted sequences):
  L2tot_D : #(i<j<k, c1, c2 in [-8,8]\\0, |d| <= D) with b_k = c1 b_i + c2 b_j + d;   D = 0, 2
  L1tot_D : #(i<k, c in 1..8, |d| <= D) with b_k = c b_i + d
  lag_ins_D, lag_oos_D : lag pattern b_k = c1 b_{k-r} + c2 b_{k-s} + d chosen on the in-sample k, hits on the rest
  mul_tot : #(i<=j, l not in {i,j}, k > all) with b_k = round(b_i b_j / b_l)
  pow_exact : #(i<k, p/q > 1, q <= 4, p <= 8q) with b_k = round(b_i^(p/q))
Output: other.json {seq: {stat: {'side': 'hi', 'real': v, 'null': [..]}}} and other_real.txt
"""
import json, sys, time, math
import numpy as np
from math import ceil

D = np.load('data.npz')
NAMES = [str(x) for x in D['names']]
PBM = (1 << 61) - 1


# ---------------- closed forms ----------------
def strips(a, fam):
    a = a.astype(np.float64)
    nn = np.arange(1, len(a) + 1.0)
    if fam == 'F1': return np.log(a), np.log(a + 1), nn
    if fam == 'F3': return np.log(a), np.log(a + 1), np.log(nn)
    if fam == 'F4': return np.log(a - .5), np.log(a + .5), nn


def best_hits(a, fam, nin):
    lo, hi, phi = strips(a, fam)
    m = len(a)
    I, J = np.triu_indices(m, 1)
    ok = phi[J] != phi[I]
    I, J = I[ok], J[ok]
    U, V = [], []
    for bi in (lo, hi):
        for bj in (lo, hi):
            v = (bj[J] - bi[I]) / (phi[J] - phi[I]); U.append(bi[I] - phi[I] * v); V.append(v)
    U = np.concatenate(U); V = np.concatenate(V)
    best, arg, oi, oo = 0, None, -1, 0
    for du in (-1, 1):
        for dv in (-1, 1):
            uu = U + du * 1e-13 * np.maximum(1, np.abs(U)); vv = V + dv * 1e-13 * np.maximum(1, np.abs(V)) / 60
            val = uu[:, None] + phi[None, :] * vv[:, None]
            ins = (val >= lo[None, :]) & (val < hi[None, :])
            c = ins.sum(1); ci = ins[:, :nin].sum(1); co = ins[:, nin:].sum(1)
            mi = int(ci.max()); mo = int(co[ci == mi].max())
            if mi > oi or (mi == oi and mo > oo): oi, oo = mi, mo
            k = int(c.argmax())
            if c[k] > best: best, arg = int(c[k]), (uu[k], vv[k])
    return best, arg, oi, oo


def exact_check(a, fam, uv):
    u, v = uv
    nn = np.arange(1, len(a) + 1, dtype=np.float64)
    if fam == 'F1': g = np.floor(np.exp(u + nn * v))
    if fam == 'F3': g = np.floor(np.exp(u) * nn ** v)
    if fam == 'F4': g = np.rint(np.exp(u + nn * v))
    return [int(i + 1) for i in np.nonzero(g == a)[0]]


# ---------------- recurrences ----------------
def bm_modp(s, p=PBM):
    s = [int(x) % p for x in s]
    C, Bc = [1], [1]
    L, m, b = 0, 1, 1
    Ls, zero = [], []
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            d = (d + C[i] * s[n - i]) % p
        zero.append(d == 0)
        if d == 0:
            m += 1
        else:
            coef = d * pow(b, p - 2, p) % p
            T = C[:]
            if len(C) < len(Bc) + m: C = C + [0] * (len(Bc) + m - len(C))
            for i in range(len(Bc)): C[i + m] = (C[i + m] - coef * Bc[i]) % p
            if 2 * L <= n:
                L = n + 1 - L; Bc = T; b = d; m = 1
            else:
                m += 1
        Ls.append(L)
    return Ls, zero


def bm_stats(a):
    Ls, zero = bm_modp(a)
    return sum(zero[m] for m in range(8, len(a))), max(ceil(m / 2) - Ls[m - 1] for m in range(8, len(a) + 1))


def grid(k, cmax=10):
    r = np.arange(-cmax, cmax + 1, dtype=np.int64)
    return np.array(np.meshgrid(*([r] * k), indexing='ij')).reshape(k, -1).T


GRIDS = {k: grid(k) for k in range(1, 4)}


def runs(M):
    rows, T = M.shape
    cur = np.zeros(rows, np.int32); best = np.zeros(rows, np.int32); ending = np.zeros((rows, T), np.int32)
    for t in range(T):
        cur = np.where(M[:, t], cur + 1, 0); best = np.maximum(best, cur); ending[:, t] = cur
    lead = np.argmin(np.concatenate([M, np.zeros((rows, 1), bool)], axis=1), axis=1)
    return best, lead, ending


def smallint_stats(a, n_in):
    a = np.array(a, np.int64)
    res = {}
    for const in (False, True):
        best_any = best_lead = 0
        choice = None
        for k in range(1, 4):
            G = GRIDS[k]
            T = len(a) - k
            lag = np.array([a[k - 1 - i: k - 1 - i + T] for i in range(k)])
            target = a[k:]
            t_end = n_in - 1 - k
            if t_end < 1: continue
            R = target[None, :] - G @ lag
            if const:
                Mq = R[:, 1:] == R[:, :-1]
                bst, lead, ending = runs(Mq)
                best_any = max(best_any, int(bst.max()) + 1); best_lead = max(best_lead, int(lead.max()) + 1)
                ins = ending[:, t_end - 1] + 1
                j = int(np.argmax(ins)); cand = int(ins[j]); dd = R[j, t_end]
                oos = int(np.sum(R[j, t_end + 1:] == dd))
            else:
                Mq = R == 0
                bst, lead, ending = runs(Mq)
                best_any = max(best_any, int(bst.max())); best_lead = max(best_lead, int(lead.max()))
                ins = ending[:, t_end]
                j = int(np.argmax(ins)); cand = int(ins[j])
                oos = int(np.sum(R[j, t_end + 1:] == 0))
            key = (cand - k, -k)
            if choice is None or key > choice[0]:
                choice = (key, k, tuple(int(x) for x in G[j]), cand, oos)
        tag = 'c' if const else 'h'
        res['si_run_' + tag] = best_any; res['si_lead_' + tag] = best_lead
        res['si_oos_' + tag] = choice[4] if choice else 0
        res['si_choice_' + tag] = choice
    return res


# ---------------- relations ----------------
CS = np.array([c for c in range(-8, 9) if c != 0], np.int64)
C1 = np.repeat(CS, CS.size); C2 = np.tile(CS, CS.size)
DS = (0, 2)


def lin_stats(a):
    a = np.asarray(a, np.int64); n = len(a)
    out = {}
    # pairs i < j ordered by j, flattened with the coefficient pairs; key = V*BIG + position
    I, J = np.triu_indices(n, 1)
    o = np.lexsort((I, J)); I, J = I[o], J[o]
    V = (a[I][:, None] * C1[None, :] + a[J][:, None] * C2[None, :]).ravel()
    BIG = 1 << 20
    keys = np.sort(V * BIG + np.arange(V.size))
    npair_lt = np.searchsorted(J, np.arange(n), 'left') * CS.size ** 2  # #flat entries with j < k
    for Dv in DS:
        tot = 0
        for d in range(-Dv, Dv + 1):
            x = a - d
            tot += int((np.searchsorted(keys, x * BIG + npair_lt) - np.searchsorted(keys, x * BIG)).sum())
        out['L2tot_D%d' % Dv] = tot
        t1 = 0
        for c in range(1, 9):
            diff = np.abs(a[None, :] - c * a[:, None])  # [i, k]
            t1 += int(((diff <= Dv) & np.triu(np.ones((n, n), bool), 1)).sum())
        out['L1tot_D%d' % Dv] = t1
    return out


def lag_oos(a, Dv, n_in, smax=6):
    a = np.asarray(a, np.int64); n = len(a)
    best = None
    for r in range(1, smax):
        for s in range(r + 1, smax + 1):
            if s >= n: continue
            k = np.arange(s, n)
            V = a[k - r][:, None] * C1[None, :] + a[k - s][:, None] * C2[None, :]
            hit = np.abs(a[k][:, None] - V) <= Dv
            ins = hit[k < n_in].sum(0); oos = hit[k >= n_in].sum(0)
            j = int(np.argmax(ins))
            if best is None or ins[j] > best[0]:
                best = (int(ins[j]), int(oos[j]), r, s, int(C1[j]), int(C2[j]))
    return best


def mul_tot(a):
    a = np.asarray(a, np.int64); n = len(a)
    I, J = np.triu_indices(n - 1, 0)
    prods = a[I] * a[J]
    BIG2 = 64
    keys = np.sort(a * BIG2 + np.arange(n))
    tot = 0
    for l in range(n - 1):
        ok = (I != l) & (J != l)
        y = int(a[l]); v = (2 * prods[ok] + y) // (2 * y)
        m = np.maximum(np.maximum(I[ok], J[ok]), l)
        tot += int((np.searchsorted(keys, v * BIG2 + BIG2) - np.searchsorted(keys, v * BIG2 + m + 1)).sum())
    return tot


FR = sorted({(p, q) for q in range(1, 5) for p in range(q + 1, 8 * q + 1) if math.gcd(p, q) == 1})


def pow_exact(a):
    af = np.asarray(a, np.float64); n = len(a)
    up = np.triu(np.ones((n, n), bool), 1)
    ex = 0
    for p, q in FR:
        v = np.round(af ** (p / q))
        ex += int(((v[:, None] == af[None, :]) & up).sum())
    return ex


def all_stats(a, nin):
    out = {}
    for fam in ('F1', 'F3', 'F4'):
        b, _, oi, oo = best_hits(a, fam, nin)
        out[fam + '_hits'] = b; out[fam + '_oos'] = oo
    out['bm_zero'], out['bm_def'] = bm_stats(a)
    si = smallint_stats(a, nin)
    for k in ('si_run_h', 'si_lead_h', 'si_oos_h', 'si_run_c', 'si_lead_c', 'si_oos_c'):
        out[k] = si[k]
    out.update(lin_stats(a))
    for Dv in DS:
        b = lag_oos(a, Dv, nin)
        out['lag_ins_D%d' % Dv] = b[0]; out['lag_oos_D%d' % Dv] = b[1]
    out['mul_tot'] = mul_tot(a)
    out['pow_exact'] = pow_exact(a)
    return out, si


if __name__ == '__main__':
    t0 = time.time()
    res = {}
    rep = []
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 9
    for nm in NAMES:
        M = D[nm]; I = D[nm + '_ins']; Ln = D[nm + '_len']
        rows = []
        for r in range(min(M.shape[0], nmax + 1)):
            a = M[r, :Ln[r]]; nin = int(I[r].sum())
            st, si = all_stats(a, nin)
            rows.append(st)
            if r == 0:
                rep.append(f'== {nm} (real, {len(a)} terms, {nin} in-sample)')
                rep.append('   ' + ', '.join(f'{k}={v}' for k, v in st.items()))
                for fam in ('F1', 'F3', 'F4'):
                    b, uv, oi, oo = best_hits(a, fam, nin)
                    par = f'u = {uv[0]:.10g}, v = {uv[1]:.10g}'
                    rep.append(f'   {fam}: {b} hits at n = {exact_check(a, fam, uv)} ({par}); OOS: in-sample best {oi}, max OOS over ties {oo}')
                rep.append(f'   small-int recurrence choice h={si["si_choice_h"]} c={si["si_choice_c"]}')
                for Dv in DS:
                    rep.append(f'   lag pattern D={Dv}: (ins, oos, r, s, c1, c2) = {lag_oos(a, Dv, nin)}')
        res[nm] = {k: {'side': 'hi', 'real': rows[0][k], 'null': [x[k] for x in rows[1:]]} for k in rows[0]}
        print(nm, 'done %.1f s' % (time.time() - t0), flush=True)
    json.dump(res, open('other.json', 'w'))
    open('other_real.txt', 'w').write('\n'.join(rep) + '\n')
    print('\n'.join(rep))
