"""Experiment 3: the Ulam spiral (1 at the centre, unit steps).

Statistics per set (larger = more striking):
  T, Q, Mx : collinear triples, 4-sets, most points on one line (any slope)
  A0       : points on the four lines through the centre (x=0, y=0, x=y, x=-y)
  Dmax     : most points on one diagonal line (x-y=c or x+y=c)
  Dpairs   : pairs of points on a common diagonal line
  Amax     : most points on one horizontal or vertical line
  Apairs   : pairs of points on a common horizontal or vertical line
Exchangeable rank p over the N+1 sets, family-wise min p and max z over the
eight statistics.   Usage: spiral.py <view> <Nnull>
"""
import sys, os, json
import numpy as np
from common import EXPS, HERE


def ulam_xy(n):
    n = np.asarray(n, dtype=np.int64)
    s = np.floor(np.sqrt((n - 1).astype(np.float64))).astype(np.int64)
    s[(s + 1) ** 2 <= n - 1] += 1
    s[s * s > n - 1] -= 1                    # s = isqrt(n-1)
    k = (s + 1) // 2
    t = 2 * k + 1
    m = t * t
    t = t - 1
    x = np.empty_like(n); y = np.empty_like(n)
    c1 = n >= m - t
    x[c1] = (k - (m - n))[c1]; y[c1] = -k[c1]
    m2 = m - t
    c2 = ~c1 & (n >= m2 - t)
    x[c2] = -k[c2]; y[c2] = (-k + (m2 - n))[c2]
    m3 = m2 - t
    c3 = ~c1 & ~c2 & (n >= m3 - t)
    x[c3] = (-k + (m3 - n))[c3]; y[c3] = k[c3]
    c4 = ~c1 & ~c2 & ~c3
    x[c4] = k[c4]; y[c4] = (k - (m3 - n - t))[c4]
    return x, y


def selftest():
    n = np.arange(1, 200001)
    x, y = ulam_xy(n)
    assert x[0] == 0 and y[0] == 0
    step = np.abs(np.diff(x)) + np.abs(np.diff(y))
    assert (step == 1).all()
    assert len(set(zip(x.tolist(), y.tolist()))) == len(n)
    # ring k holds (2k-1)^2+1 .. (2k+1)^2, all with max(|x|,|y|) = k
    k = (np.floor(np.sqrt(n - 1)) + 1) // 2
    assert (np.maximum(np.abs(x), np.abs(y)) == k).all()


def pair_count_max(keys_by_set):
    """keys_by_set: (S, K) int64 line labels; returns (max per set, pairs per set)."""
    S, K = keys_by_set.shape
    srt = np.sort(keys_by_set, axis=1)
    st = np.ones_like(srt, dtype=bool)
    st[:, 1:] = srt[:, 1:] != srt[:, :-1]
    runid = np.cumsum(st.ravel()) - 1
    cnt = np.bincount(runid)
    D = st.sum(axis=1)
    off = np.concatenate([[0], np.cumsum(D)[:-1]])
    mx = np.maximum.reduceat(cnt, off)
    pairs = np.add.reduceat(cnt * (cnt - 1) // 2, off)
    return mx, pairs


def stats(X):
    S, K = X.shape
    x, y = ulam_xy(X.ravel())
    x = x.reshape(S, K); y = y.reshape(S, K)
    A0 = ((x == 0) | (y == 0) | (x == y) | (x == -y)).sum(axis=1)
    BIG = 1 << 20
    dmx1, dp1 = pair_count_max(x - y)
    dmx2, dp2 = pair_count_max(x + y)
    amx1, ap1 = pair_count_max(x)
    amx2, ap2 = pair_count_max(y)
    Dmax = np.maximum(dmx1, dmx2); Dpairs = dp1 + dp2
    Amax = np.maximum(amx1, amx2); Apairs = ap1 + ap2
    # collinear, any slope, in chunks of sets
    I, J = np.triu_indices(K, 1)
    T = np.zeros(S); Q = np.zeros(S); Mx = np.zeros(S, dtype=np.int64)
    R = int(max(np.abs(x).max(), np.abs(y).max()))
    SB = int(np.ceil(np.log2(4 * R + 2)))
    for c0 in range(0, S, 500):
        xs, ys = x[c0:c0 + 500], y[c0:c0 + 500]
        s_ = xs.shape[0]
        dx = (xs[:, J] - xs[:, I]).ravel()
        dy = (ys[:, J] - ys[:, I]).ravel()
        g = np.gcd(dx, dy)
        a = dx // g; b = dy // g
        neg = (a < 0) | ((a == 0) & (b < 0))
        a[neg] = -a[neg]; b[neg] = -b[neg]
        SI = (np.arange(s_)[:, None] * K + I[None, :]).ravel()
        key = (SI << (2 * SB)) | (a << SB) | (b + 2 * R)
        key.sort()
        st = np.empty(len(key), dtype=bool); st[0] = True
        st[1:] = key[1:] != key[:-1]
        cnt = np.bincount(np.cumsum(st) - 1).astype(np.float64)
        gset = (key[st] >> (2 * SB)) // K
        T[c0:c0 + s_] = np.bincount(gset, cnt * (cnt - 1) / 2, minlength=s_)
        Q[c0:c0 + s_] = np.bincount(gset, cnt * (cnt - 1) * (cnt - 2) / 6, minlength=s_)
        off = np.searchsorted(gset, np.arange(s_), side="left")
        Mx[c0:c0 + s_] = np.maximum.reduceat(cnt, off).astype(np.int64) + 1
    return {"T": T, "Q": Q, "Mx": Mx, "A0": A0, "Dmax": Dmax, "Dpairs": Dpairs,
            "Amax": Amax, "Apairs": Apairs}


if __name__ == "__main__":
    selftest()
    view = sys.argv[1]
    NN = int(sys.argv[2])
    null = np.load(os.path.join(HERE, f"null_{view}.npy"))[:NN]
    lim = {"full": 10**9, "page": 500000, "h5": 100000}[view]
    real = np.array([e for e in EXPS if e < lim], dtype=np.int64)
    X = np.vstack([real[None, :], null])
    S = X.shape[0]
    st = stats(X)
    out = {"view": view, "N": S - 1, "stats": {}}
    minp = np.ones(S); maxz = np.full(S, -np.inf)
    for name, s in st.items():
        s = np.asarray(s, dtype=np.float64)
        sa = np.sort(s)
        pall = (S - np.searchsorted(sa, s, side="left")) / S
        mu, sd = s[1:].mean(), s[1:].std()
        z = (s - mu) / sd if sd > 0 else np.zeros_like(s)
        minp = np.minimum(minp, pall); maxz = np.maximum(maxz, z)
        out["stats"][name] = {"real": float(s[0]), "null_mean": float(mu),
                              "null_sd": float(sd), "p": float(pall[0]),
                              "z": float(z[0])}
    out["familywise_p_minp"] = int((minp <= minp[0] + 1e-12).sum()) / S
    out["familywise_p_maxz"] = int((maxz >= maxz[0] - 1e-12).sum()) / S
    x, y = ulam_xy(real)
    out["coords"] = {int(e): (int(a), int(b)) for e, a, b in zip(real, x, y)}
    with open(os.path.join(HERE, f"spiral_{view}.json"), "w") as f:
        json.dump(out, f, indent=1)
    for k, v in out["stats"].items():
        print(k, {a: round(b, 4) for a, b in v.items()})
    print("familywise min-p:", out["familywise_p_minp"], " max-z:", out["familywise_p_maxz"])
