#!/usr/bin/env python3
"""Mersenne population for the transfer: primes 2 <= k < B = 83,880,077 (b = 2, kind 'm').

Sample: every prime k < 2^16 (weight 1), and in each octave [2^i, 2^(i+1)) ∩ [.., B), i = 16..26,
NS = 3000 primes drawn uniformly without replacement among the non-exponents (weight = count/NS).
The 51 exponents (weight 1) and the 52nd, 136,279,841 (scored only), are added.
Trial factoring to depth K = 10^4: s(k) = 1 if no q = 2jk + 1, j <= K, q = +-1 mod 8, q < 2^k - 1,
divides 2^k - 1 (vectorised; 2^k mod q with an exact float-quotient mulmod, q < 2^42).
Features: feats.features('m', 2, k), the training thresholds and standardisation.
Output: mersenne.npz."""
import json, math, os, sys, time
import numpy as np
import feats as F

HERE = os.path.dirname(os.path.abspath(__file__))
B = 83880077
K = 10000
NS = 3000
SEED = 52
EXP = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
       9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
       1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951,
       30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
NEW = 136279841
SMALL = [3, 5, 7, 11, 13, 17, 19, 23]


def odd_primes_below(n):
    s = np.ones(n // 2, bool); s[0] = False
    for i in range(3, int(n ** 0.5) + 1, 2):
        if s[i // 2]:
            s[i * i // 2::i] = False
    return 2 * np.flatnonzero(s).astype(np.int64) + 1


def mulmod(a, b, q, qf):
    quo = np.floor(a.astype(np.float64) * b.astype(np.float64) / qf).astype(np.uint64)
    r = (a * b - quo * q).view(np.int64)
    qi = q.view(np.int64)
    r = np.where(r < 0, r + qi, r)
    r = np.where(r >= qi, r - qi, r)
    return r.view(np.uint64)


def pow2mod(e, q):
    """2^e mod q elementwise (e, q uint64; q < 2^42)."""
    qf = q.astype(np.float64)
    r = np.ones_like(q)
    nb = int(e.max()).bit_length()
    for i in range(nb - 1, -1, -1):
        r = mulmod(r, r, q, qf)
        bit = ((e >> np.uint64(i)) & np.uint64(1)).astype(bool)
        r2 = r << np.uint64(1)
        r2 = np.where(r2 >= q, r2 - q, r2)
        r = np.where(bit, r2, r)
    return r


def tf_survive(ks, K=K, batch=200):
    """s(k) for each k (array of primes)."""
    surv = np.ones(len(ks), bool)
    j = np.arange(1, K + 1, dtype=np.uint64)
    for st in range(0, len(ks), batch):
        kk = ks[st:st + batch].astype(np.uint64)
        q = (np.uint64(2) * kk[:, None] * j[None, :] + np.uint64(1))
        e = np.broadcast_to(kk[:, None], q.shape)
        row = np.broadcast_to(np.arange(len(kk))[:, None], q.shape)
        m8 = q % np.uint64(8)
        keep = (m8 == 1) | (m8 == 7)
        for p in SMALL:
            keep &= (q % np.uint64(p) != 0) | (q == np.uint64(p))
        # q must be a proper divisor of 2^k - 1
        small = kk < 64
        if small.any():
            lim = np.array([(1 << int(x)) - 1 if x < 64 else 2 ** 64 - 1 for x in kk], dtype=np.uint64)
            keep &= q < lim[:, None]
        qv, ev, rv = q[keep], e[keep], row[keep]
        if len(qv) == 0:
            continue
        r = pow2mod(ev, qv)
        hit = np.unique(rv[r == 1])
        surv[st + hit] = False
    return surv


def main():
    t0 = time.time()
    P = odd_primes_below(B)
    P = np.concatenate(([2], P))
    expset = set(EXP)
    isexp = np.isin(P, EXP)
    assert isexp.sum() == 51
    rng = np.random.default_rng(SEED)
    ks, W, stratum = [], [], []
    lo = P[(P < 2 ** 16) & ~isexp]
    ks.append(lo); W.append(np.ones(len(lo))); stratum.append(np.full(len(lo), 15))
    for i in range(16, 27):
        a, b = 2 ** i, min(2 ** (i + 1), B)
        pool = P[(P >= a) & (P < b) & ~isexp]
        pick = np.sort(rng.choice(pool, NS, replace=False))
        ks.append(pick); W.append(np.full(NS, len(pool) / NS)); stratum.append(np.full(NS, i))
    ks.append(np.array(EXP, np.int64)); W.append(np.ones(51)); stratum.append(np.full(51, -1))
    ks.append(np.array([NEW], np.int64)); W.append(np.zeros(1)); stratum.append(np.full(1, -2))
    ks = np.concatenate(ks).astype(np.int64); W = np.concatenate(W); stratum = np.concatenate(stratum)
    print(f"population {len(P)} primes, sample {len(ks)} rows, weight sum {W.sum():.0f}", flush=True)
    del P
    surv = tf_survive(ks)
    print(f"TF done {time.time() - t0:.0f}s; survivors {surv.sum()} of {len(ks)}; "
          f"exponents all survive: {surv[stratum < 0].all()}", flush=True)
    # cross-check with the earlier kmin table (primes < 10^5, depth 10^4) if present
    km = os.path.join(HERE, "..", "kmin_100000_10000.npy")
    if os.path.exists(km):
        A = np.load(km)
        ref = dict(zip(A[0].tolist(), (A[1] == 0).tolist()))
        m = [(int(k), bool(s)) for k, s in zip(ks, surv) if int(k) in ref]
        bad = [k for k, s in m if ref[k] != s]
        print(f"cross-check vs kmin_100000_10000.npy: {len(m)} primes, {len(bad)} disagreements {bad[:10]}", flush=True)
    d = np.load(os.path.join(HERE, "train.npz"))
    names_train = list(d["names"]); thr = json.loads(str(d["thr"]))
    # only survivors carry weight (w = 0 otherwise); features in chunks to bound memory
    sel = surv | (stratum < 0)
    ks, W, stratum, surv = ks[sel], W[sel], stratum[sel], surv[sel]
    parts = []
    for c0 in range(0, len(ks), 4000):
        X, names, groups, types = F.features("m", 2, ks[c0:c0 + 4000])
        X, names, groups, types, _ = F.add_thresholds(X, names, groups, types, thr)
        idx = [names.index(nm) for nm in names_train]
        parts.append(((X[:, idx] - d["mu"]) / d["sd"]).astype(np.float32))
        del X
    X = np.vstack(parts)
    np.savez_compressed(os.path.join(HERE, "mersenne.npz"), X=X.astype(np.float32), k=ks, W=W,
                        stratum=stratum, surv=surv)
    print(f"features {X.shape} {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
