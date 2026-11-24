"""The class-response law (derivation.md section 6) evaluated for round 83 (bands of D_k) and rounds 85-86 (last-digit
tables), retrospective.

Short-range model: the sieve model M_Q, Q = 30030: every residue coprime to Q is 'prime' independently with probability
p = (Q/phi(Q))/lam. In M_Q the class counts N_c (c mod Q) are independent with Var = (1 - p) E N_c, and by periodicity
the linear (first-order Hoeffding) part of any window statistic is exactly its regression on the N_c. The law:

  E chi^2(primes) = E chi^2(M_Q) - FO_sieve + FO_HL,
  FO_X = sum_cells (1/E_cell) sum_{c,c'} r_{cell,c} Cov_X(N_c, N_c') r_{cell,c'},

r = response of the cell residual (count minus expected) to one more prime in class c (per site), Cov_sieve =
diag (1-p) E N_c, Cov_HL = closed form of section 4 at modulus Q and window length H.

Responses are measured by insertion: for each non-prime admissible site m of a simulated block, the change of every
cell count if m were prime (this change does not depend on x_m, so averaging over non-prime sites is unbiased).
Two halves of the blocks give two independent estimates r1, r2; FO uses r1 r2 (no noise bias).
E chi^2(M_Q) is the block mean of chi^2 against the pooled law, x B/(B-1).
usage: classresp.py lam_x window_H DMAX blocks span seed out.json"""
import sys, math, json, time
import numpy as np

X = float(sys.argv[1]); HWIN = float(sys.argv[2]); DMAX = int(sys.argv[3])
NB = int(sys.argv[4]); SPAN = int(float(sys.argv[5])); SEED = int(sys.argv[6]); OUT = sys.argv[7]
LAM = math.log(X)
Q = 30030; PQ = (2, 3, 5, 7, 11, 13)
adm = np.array([a for a in range(Q) if math.gcd(a, Q) == 1], np.int64); NC = len(adm)
cidx = np.full(Q, -1, np.int64); cidx[adm] = np.arange(NC)
p = (Q / NC) / LAM
EG = 0.5772156649015329; BMS = 1 - EG - math.log(2 * math.pi)
SQ = sum(math.log(r) / (r - 1) for r in PQ)
vQ = 1 - (math.log(HWIN / Q) + SQ - BMS) / LAM
dig = np.full(10, -1); dig[[1, 3, 7, 9]] = np.arange(4)
rng = np.random.default_rng(SEED)

# statistics: ('band', k) cells 0..19 ; ('dig', k) cells a*4+b
STATS = [('band', k) for k in (1, 2, 3, 4)] + [('dig', k) for k in (1, 2, 3)]
NCELL = {'band': 20, 'dig': 16}


def cell_of(kind, left, right):
    d = right - left
    ok = d <= (240 if kind == 'band' else DMAX)
    if kind == 'band':
        c = (d - 2) // 12
        ok &= (c >= 0) & (c < 20)
    else:
        c = dig[left % 10] * 4 + dig[right % 10]
    return np.where(ok, c, -1)


def block_counts(P, kind, k):
    c = cell_of(kind, P[:-k], P[k:])
    c = c[c >= 0]
    return np.bincount(c, minlength=NCELL[kind]).astype(float)


def insertion(P, M, kind, k, acc):
    """add to acc[class(m), cell] the change of cell counts when m is inserted (m not prime)."""
    j = np.searchsorted(P, M) - 1
    ok = (j - k + 1 >= 0) & (j + k < len(P))
    M = M[ok]; j = j[ok]; cl = cidx[M % Q]
    nc = NCELL[kind]
    def put(left, right, w):
        c = cell_of(kind, left, right)
        sel = c >= 0
        np.add(acc, w * np.bincount(cl[sel] * nc + c[sel], minlength=NC * nc).reshape(NC, nc), out=acc)
    for t in range(k):                                   # removed pairs (P[i], P[i+k]), i = j-t
        i = j - t
        put(P[i], P[i + k], -1.0)
    put(P[j - k + 1], M, 1.0); put(M, P[j + k], 1.0)      # added pairs ending / starting at m
    for i_off in range(k - 1):                           # added (P[i], P[i+k-1]), i = j-k+2 .. j
        i = j - k + 2 + i_off
        put(P[i], P[i + k - 1], 1.0)


def main():
    t0 = time.time()
    counts = {s: [] for s in STATS}
    acc = {s: [np.zeros((NC, NCELL[s[0]])), np.zeros((NC, NCELL[s[0]]))] for s in STATS}
    nsite = [np.zeros(NC), np.zeros(NC)]
    periods = 0
    for b in range(NB):
        base = (np.arange(SPAN // Q, dtype=np.int64)[:, None] * Q + adm[None, :]).ravel() + Q   # start at Q
        isp = rng.random(len(base)) < p
        P = base[isp]; M = base[~isp]
        periods += SPAN // Q
        h = b % 2
        nsite[h] += np.bincount(cidx[M % Q], minlength=NC)
        for s in STATS:
            counts[s].append(block_counts(P, *s))
            insertion(P, M, s[0], s[1], acc[s][h])
        del base, isp, P, M
        if b % 5 == 0:
            print('block', b, '%.0f s' % (time.time() - t0), flush=True)
    res = {'lam': LAM, 'x': X, 'p': p, 'H_window': HWIN, 'v_Q': vQ, 'blocks': NB, 'span': SPAN, 'stats': {}}
    # fibre structure for the cross-class covariance: classes that differ only mod r (r | Q)
    fib = {r: (adm % (Q // r)) for r in PQ}
    for s in STATS:
        kind, k = s
        C = np.array(counts[s])                              # (NB, ncell)
        pooled = C.sum(0)
        if kind == 'band':
            Pc = pooled / pooled.sum()
            E = C.sum(1, keepdims=True) * Pc[None, :]
            keepc = (Pc * 4.8e6 >= 50) if k == 1 else np.ones(20, bool)   # round 83: N = 4.8e6 per k; all 20 bands for k >= 2
        else:
            T = pooled.reshape(4, 4); Pab = (T / T.sum(1, keepdims=True)).ravel()
            rows = C.reshape(NB, 4, 4).sum(2)
            E = np.repeat(rows, 4, axis=1) * Pab[None, :]
            keepc = np.ones(16, bool)
        chi = (((C - E) ** 2 / E)[:, keepc]).sum(1)
        dshort = float(chi.mean() * NB / (NB - 1)); dsd = float(chi.std(ddof=1))
        # responses per site (two halves); residual response
        r = []
        for h in (0, 1):
            psi = acc[s][h] / nsite[h][:, None]
            if kind == 'band':
                tot = psi.sum(1, keepdims=True)              # response of the number of pairs with D <= 240
                rr = psi - tot * Pc[None, :]
            else:
                rowtot = psi.reshape(NC, 4, 4).sum(2)        # response of the row totals N_a
                rr = (psi.reshape(NC, 4, 4) - rowtot[:, :, None] * Pab.reshape(4, 4)[None]).reshape(NC, 16)
            r.append(rr[:, keepc])
        Ecell = (pooled / periods)[keepc]                     # expected cell count per period of Q (one site per class)
        FO_s = float(np.sum((p * (1 - p)) * (r[0] * r[1]).sum(0) / Ecell))
        # HL covariance per period: p [ v_Q delta - sum_r ln r/((r-1) lam) 1(c' in fibre_r(c), c' != c) ]
        cross = np.zeros_like(r[1])
        for q_ in PQ:
            fs = np.zeros((Q // q_, r[1].shape[1]))
            np.add.at(fs, fib[q_], r[1])
            cross += (math.log(q_) / ((q_ - 1) * LAM)) * (fs[fib[q_]] - r[1])
        FO_hl = float(np.sum(p * (r[0] * (vQ * r[1] - cross)).sum(0) / Ecell))
        dof = int(keepc.sum()) - (1 if kind == 'band' else 4)
        res['stats']['%s_k%d' % s] = {'dof': dof, 'chi2_sieve_model': dshort, 'chi2_sieve_sd_per_block': dsd,
                                      'FO_sieve': FO_s, 'FO_HL': FO_hl, 'chi2_primes_law': dshort - FO_s + FO_hl}
        print(s, {a: round(v, 3) if isinstance(v, float) else v for a, v in res['stats']['%s_k%d' % s].items()}, flush=True)
    json.dump(res, open(OUT, 'w'), indent=1)
    print('done %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
