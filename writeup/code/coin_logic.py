import numpy as np, lzma, bz2, zlib, subprocess, sys
from math import log, log2, sqrt

X = 10**9
W = 10**8
GAMMA = 0.5772156649015329


def small_primes(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def h2(p):
    return 0.0 if p <= 0 or p >= 1 else -(p * log2(p) + (1 - p) * log2(1 - p))


def cond_entropy(keys, vals):
    order = np.lexsort((vals, keys))
    k, v = keys[order], vals[order]
    kv = k.astype(np.int64) * (vals.max() + 1) + v
    _, cnt_kv = np.unique(kv, return_counts=True)
    _, cnt_k = np.unique(k, return_counts=True)
    n = len(keys)
    return (-(cnt_kv / n * np.log2(cnt_kv / n)).sum()) - (-(cnt_k / n * np.log2(cnt_k / n)).sum())


def hit_rate(keys_tr, vals_tr, keys_te, vals_te):
    best = {}
    kv = keys_tr.astype(np.int64) * (vals_tr.max() + 1) + vals_tr
    u, c = np.unique(kv, return_counts=True)
    for key, val, cnt in zip(u // (vals_tr.max() + 1), u % (vals_tr.max() + 1), c):
        if cnt > best.get(key, (0, -1))[0]:
            best[key] = (cnt, val)
    default = np.bincount(vals_tr).argmax()
    pred = np.array([best.get(k, (0, default))[1] for k in keys_te])
    return (pred == vals_te).mean()


top = int(sqrt(X + W + 2)) + 1
primes = small_primes(top)
comp = np.zeros(W + 2, dtype=bool)
for p in primes:
    comp[(-X) % p::p] = True
P = ~comp
y = int(X**0.405)
surv = np.zeros(W + 2, dtype=bool)
surv[:] = True
for p in primes[primes <= y]:
    surv[(-X) % p::p] = False
idx = np.flatnonzero(surv[:W])
coin = P[idx]
rho = coin.mean()
nbits = coin.size
print(f"window [{X:.0e}, {X:.0e}+{W:.0e}], sieve to y = x^0.405 = {y}: {nbits} survivors, coin rho = {rho:.4f}, entropy bound h(rho) = {h2(rho):.4f} bits/survivor = {h2(rho)*nbits/8/1024:.0f} KB")
packed = np.packbits(coin)
raw = packed.tobytes()
print("compressibility of the coin sequence (bytes; the bound is the entropy):")
print(f"  raw packed bits        {len(raw):9d}")
print(f"  entropy bound          {int(h2(rho)*nbits/8):9d}")
print(f"  zlib -9                {len(zlib.compress(raw, 9)):9d}")
print(f"  bz2 -9                 {len(bz2.compress(raw, 9)):9d}")
print(f"  xz -9e                 {len(lzma.compress(raw, preset=9 | lzma.PRESET_EXTREME)):9d}")
rng = np.random.default_rng(1)
shuf = np.packbits(rng.permutation(coin)).tobytes()
print(f"  xz -9e, shuffled coin  {len(lzma.compress(shuf, preset=9 | lzma.PRESET_EXTREME)):9d}")
try:
    out = subprocess.run(["zstd", "-q", "--ultra", "-22", "-c"], input=raw, capture_output=True, timeout=100).stdout
    print(f"  zstd --ultra -22       {len(out):9d}")
except Exception as e:
    print("  zstd not available:", e)

fails = idx[~coin] + X
print()
print("periods hidden in the coin's failures (composite survivors): a modulus q in (y, 2y] divides some failure iff q is prime")
qs = np.arange(y + 1, 2 * y + 1)
hits = np.array([(fails % q == 0).any() for q in qs])
isprime_q = np.isin(qs, primes)
print(f"  moduli tested {len(qs)}, with a failure divisible by q: {hits.sum()}, primes among them: {(hits & isprime_q).sum()}, primes missed: {(isprime_q & ~hits).sum()}, composites hit: {(hits & ~isprime_q).sum()}")

print()
print("twin gaps below 1e9: what predicts the distance to the next twin pair?")
s = np.ones(X + 3, dtype=bool)
s[:2] = False
for i in range(2, int((X + 2) ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
tw = np.flatnonzero(s[:-2] & s[2:])
tw = tw[tw > 5]
g = np.diff(tw)
p = tw[:-1]
cap = 6 * 40
gb = np.minimum(g, cap) // 6
half = len(g) // 2
feats = {
    "none": np.zeros(len(g), dtype=np.int64),
    "p mod 5": p % 5,
    "p mod 35": p % 35,
    "p mod 385": p % 385,
    "p mod 5005": p % 5005,
    "last gap": np.concatenate(([0], gb[:-1])),
    "last 2 gaps": np.concatenate(([0, 0], gb[:-2])) * 64 + np.concatenate(([0], gb[:-1])),
    "p mod 35 + last gap": (p % 35) * 64 + np.concatenate(([0], gb[:-1])),
    "p mod 385 + last gap": (p % 385) * 64 + np.concatenate(([0], gb[:-1])),
}
H0 = cond_entropy(feats["none"], gb)
print(f"  {len(g)} twin gaps, entropy of the gap (in units of 6, capped at {cap}): {H0:.3f} bits")
print("  feature                 H(gap | feature)   exact-hit rate (train first half, test second half)")
for name, f in feats.items():
    H = cond_entropy(f, gb)
    hr = hit_rate(f[:half], gb[:half], f[half:], gb[half:])
    print(f"  {name:22s}  {H:7.3f}            {100*hr:5.1f} %")
