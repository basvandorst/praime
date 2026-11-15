"""Angle 1: linear recurrences.

1a  Berlekamp-Massey over F_P (P = 2^61-1, equal to the rational linear complexity
    except for P dividing a Hankel minor: probability ~1e-17 per minor).
    stats: bm_zero  = number of m in 8..51 where the shortest recurrence for the first
                      m terms predicts a_{m+1} exactly (zero discrepancy)
           bm_def   = max_m (ceil(m/2) - L(m)), m = 8..52 (generic value 0)
    real only: exact BM over Q for m = 8..30 (order, coefficient height, prediction).
1b  P-recursive: sum_{i=0..k} P_i(n) a_{n+i} = 0, deg P_i <= d, k = 1..6, d = 0..4,
    prefix rows; rank tracked mod q = 2^26-5 (a 26-bit prime), deficits rechecked
    mod a second prime for the real sequence.
    stat: prec_exc = max over (k,d) of the number of dependent rows before the kernel
          dies (= number of later terms predicted by a recurrence fitted on fewer terms)
          prec_sum = sum of these over (k,d).
1c  Small integer recurrences a_n = sum_{i=1..k} c_i a_{n-i} (+ d), k <= 4, |c_i| <= 10,
    d any integer; exhaustive.
    stats: longest run of consecutive n satisfied (any position / from the start),
           out of sample: recurrence chosen on n <= 25 (longest run ending at n = 25),
           number of n in 26..52 it satisfies.
"""
import json, sys, time
from fractions import Fraction
from math import ceil
import itertools
import numpy as np

P = (1 << 61) - 1


def bm_modp(s, p=P):
    s = [x % p for x in s]
    C, B = [1], [1]
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
            if len(C) < len(B) + m:
                C = C + [0] * (len(B) + m - len(C))
            for i in range(len(B)):
                C[i + m] = (C[i + m] - coef * B[i]) % p
            if 2 * L <= n:
                L = n + 1 - L
                B = T
                b = d
                m = 1
            else:
                m += 1
        Ls.append(L)
    return Ls, zero  # Ls[m-1] = L(m); zero[n] = discrepancy zero at term n (0-based)


def bm_stats(a):
    Ls, zero = bm_modp(a)
    bm_zero = sum(zero[m] for m in range(8, len(a)))  # predicts a_{m+1} from first m
    bm_def = max(ceil(m / 2) - Ls[m - 1] for m in range(8, len(a) + 1))
    return bm_zero, bm_def


def bm_exact(s):
    """BM over Q; returns for each m the connection polynomial C (C[0]=1) and L."""
    s = [Fraction(x) for x in s]
    C, B = [Fraction(1)], [Fraction(1)]
    L, m, b = 0, 1, Fraction(1)
    out = []
    for n in range(len(s)):
        d = s[n] + sum(C[i] * s[n - i] for i in range(1, L + 1))
        if d == 0:
            m += 1
        else:
            coef = d / b
            T = C[:]
            if len(C) < len(B) + m:
                C = C + [Fraction(0)] * (len(B) + m - len(C))
            for i in range(len(B)):
                C[i + m] -= coef * B[i]
            if 2 * L <= n:
                L = n + 1 - L
                B = T
                b = d
                m = 1
            else:
                m += 1
        out.append((L, C[:L + 1] + [Fraction(0)] * max(0, L + 1 - len(C))))
    return out


Q1 = 67108859  # prime < 2^26
Q2 = 67108837  # another prime < 2^26


def prec_excess(a, k, d, q=Q1):
    a = np.array(a, dtype=object)
    n_rows = len(a) - k
    U = (k + 1) * (d + 1)
    Bm = np.zeros((0, U), dtype=np.int64)
    piv = []
    r = 0
    dep = 0
    for t in range(n_rows):
        nn = t + 1
        row = np.array([(pow(nn, e, q) * (int(a[t + i]) % q)) % q for i in range(k + 1) for e in range(d + 1)],
                       dtype=np.int64)
        if r:
            coef = row[piv]
            row = (row - (coef @ Bm) % q) % q
        nz = np.nonzero(row)[0]
        if nz.size == 0:
            dep += 1  # dependent row while kernel nontrivial (r < U always here)
            continue
        c = int(nz[0])
        inv = pow(int(row[c]), q - 2, q)
        row = (row * inv) % q
        if r:
            Bm = (Bm - np.outer(Bm[:, c], row) % q) % q
        Bm = np.vstack([Bm, row])
        piv.append(c)
        r += 1
        if r == U:
            break
    return dep


def prec_stats(a):
    ex = {}
    for k in range(1, 7):
        for d in range(0, 5):
            ex[(k, d)] = prec_excess(a, k, d)
    return max(ex.values()), sum(ex.values()), ex


def grid(k, cmax=10):
    rng = np.arange(-cmax, cmax + 1, dtype=np.int64)
    g = np.array(np.meshgrid(*([rng] * k), indexing='ij')).reshape(k, -1).T
    return g  # (21^k, k)


GRIDS = {k: grid(k) for k in range(1, 5)}


def runs(M):
    """M bool (rows, T): longest run of True per row, leading run, run ending at each col."""
    rows, T = M.shape
    cur = np.zeros(rows, dtype=np.int32)
    best = np.zeros(rows, dtype=np.int32)
    ending = np.zeros((rows, T), dtype=np.int32)
    for t in range(T):
        cur = np.where(M[:, t], cur + 1, 0)
        best = np.maximum(best, cur)
        ending[:, t] = cur
    lead = np.argmin(np.concatenate([M, np.zeros((rows, 1), bool)], axis=1), axis=1)
    return best, lead, ending


def smallint_stats(a, n_in=25):
    a = np.array(a, dtype=np.int64)
    res = {}
    for const in (False, True):
        best_any = 0
        best_lead = 0
        choice = None  # (in-sample run, -k, k, row, oos)
        for k in range(1, 5):
            G = GRIDS[k]
            T = len(a) - k
            lag = np.array([a[k - 1 - i: k - 1 - i + T] for i in range(k)])  # lag[i,t] = a_{n-1-i}
            target = a[k:]
            # equation t corresponds to (0-based) term index n = t + k; 1-based n+1
            t_end = n_in - 1 - k  # last equation with 1-based index <= n_in
            for c0 in range(0, G.shape[0], 40000):
                g = G[c0:c0 + 40000]
                R = target[None, :] - g @ lag
                if const:
                    M = R[:, 1:] == R[:, :-1]
                    best, lead, ending = runs(M)
                    eq_best = best + 1
                    eq_lead = lead + 1
                    ins = ending[:, t_end - 1] + 1  # equations ending at t_end
                    best_any = max(best_any, int(eq_best.max()))
                    best_lead = max(best_lead, int(eq_lead.max()))
                    j = int(np.argmax(ins))
                    cand = int(ins[j])
                    dd = R[j, t_end]
                    oos = int(np.sum(R[j, t_end + 1:] == dd))
                else:
                    M = R == 0
                    best, lead, ending = runs(M)
                    best_any = max(best_any, int(best.max()))
                    best_lead = max(best_lead, int(lead.max()))
                    ins = ending[:, t_end]
                    j = int(np.argmax(ins))
                    cand = int(ins[j])
                    oos = int(np.sum(R[j, t_end + 1:] == 0))
                key = (cand - k, -k)  # in-sample equations beyond the k coefficients
                if choice is None or key > choice[0]:
                    choice = (key, k, tuple(int(x) for x in g[j]), cand, oos)
        tag = 'c' if const else 'h'
        res['si_run_' + tag] = best_any
        res['si_lead_' + tag] = best_lead
        res['si_oos_' + tag] = choice[4]
        res['si_choice_' + tag] = choice
    return res


def all_stats(a):
    out = {}
    out['bm_zero'], out['bm_def'] = bm_stats(a)
    mx, sm, _ = prec_stats(a)
    out['prec_exc'] = mx
    out['prec_sum'] = sm
    si = smallint_stats(a)
    for k in ('si_run_h', 'si_lead_h', 'si_oos_h', 'si_run_c', 'si_lead_c', 'si_oos_c'):
        out[k] = si[k]
    return out, si


SIDES = {'bm_zero': 'hi', 'bm_def': 'hi', 'prec_exc': 'hi', 'prec_sum': 'hi',
         'si_run_h': 'hi', 'si_lead_h': 'hi', 'si_oos_h': 'hi',
         'si_run_c': 'hi', 'si_lead_c': 'hi', 'si_oos_c': 'hi'}

if __name__ == '__main__':
    D = json.load(open('nulls.json'))
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    t0 = time.time()
    real, si = all_stats(D['real'])
    rep = open('angle1_real.txt', 'w')
    rep.write('real stats: %s\n' % real)
    rep.write('small-int choice (in-sample n<=25): h=%s c=%s\n' % (si['si_choice_h'], si['si_choice_c']))
    # exact BM detail for the real sequence
    ex = bm_exact(D['real'])
    rep.write('\nExact BM over Q, first m terms: m, L(m), max coefficient height (bits), '
              'prediction of a_{m+1}, true, correct, run of correct predictions\n')
    A = D['real']
    for m in range(8, 31):
        L, C = ex[m - 1]
        h = max(max(abs(c.numerator).bit_length(), c.denominator.bit_length()) for c in C)
        pred = -sum(C[i] * A[m - i] for i in range(1, L + 1))
        run = 0
        for mm in range(m, len(A)):
            pr = -sum(C[i] * A[mm - i] for i in range(1, L + 1))
            if pr == A[mm]:
                run += 1
            else:
                break
        predf = float(pred)
        rep.write('%2d L=%2d height=%5d bits pred=%.6g true=%d correct=%s run=%d\n'
                  % (m, L, h, predf, A[m], pred == A[m], run))
    _, _, exd = prec_stats(D['real'])
    rep.write('\nP-recursive dependent rows per (k,d) mod Q1: %s\n' % {f'{k}': v for k, v in exd.items() if v})
    # recheck nonzero ones with Q2
    for kd, v in exd.items():
        if v:
            rep.write('recheck %s mod Q2: %d\n' % (kd, prec_excess(D['real'], kd[0], kd[1], Q2)))
    rep.close()
    out = {s: {'side': SIDES[s], 'real': real[s], 'a': [], 'b': []} for s in SIDES}
    for fam in ('a', 'b'):
        for seq in D[fam][:nmax]:
            st, _ = all_stats(seq)
            for s in SIDES:
                out[s][fam].append(st[s])
    json.dump(out, open('angle1.json', 'w'))
    print('done', time.time() - t0)
