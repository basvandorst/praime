"""Independent replay of the prime gap stacks game (same rules, same LCG) for a few settings."""
import json, sys
import numpy as np
N = 20_000_000
s = np.ones(N + 1, dtype=bool); s[:2] = False
for p in range(2, int(N ** 0.5) + 1):
    if s[p]: s[p * p::p] = False
P = np.nonzero(s)[0].astype(np.int64); del s
G = np.diff(P)
def play(W, H, rule, limit=None):
    h = [0] * W; lines = 0; seed = 20261121; colUse = [0] * W
    for i in range(len(G) if limit is None else min(limit, len(G))):
        g = int(G[i]); L = max(1, g >> 1)
        if rule == "mod": c = int(P[i]) % W
        elif rule == "rand":
            seed = (seed * 1664525 + 1013904223) % 2**32; c = int((seed / 2**32) * W)
        else: c = min(range(W), key=lambda k: (h[k], k))
        if h[c] + L > H: return {"W": W, "H": H, "rule": rule, "placed": i, "lines": lines, "over_gap": g, "over_after": int(P[i]), "col": c}
        h[c] += L; colUse[c] += 1
        m = min(h)
        if m: h = [x - m for x in h]; lines += m
    return {"W": W, "H": H, "rule": rule, "placed": i + 1, "lines": lines, "end": True}
out = [play(12, 40, "rand"), play(12, 40, "mod"), play(30, 60, "rand"), play(7, 40, "mod"), play(12, 40, "low", 300000)]
json.dump(out, open("ref40game.json", "w")); print("\n".join(map(str, out)))
