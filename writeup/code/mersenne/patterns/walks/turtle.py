"""Turtle layout: page check, runs = prime gaps, displacement, loops, heading classes,
against block-shuffled gap sequences (same gaps, order randomised inside blocks of 500)."""
import json, math
import numpy as np
from common import series, load_pos, N

S = series(N)
X, Y, head, primes = S["X"], S["Y"], S["head"], S["primes"]
out = {}

# page check
fx, fy = load_pos("turtle")
x0, x1, y0, y1 = X[1:].min(), X[1:].max(), Y[1:].min(), Y[1:].max()
sc = 1400 / max(x1 - x0, y1 - y0, 1)
out["page_check"] = dict(max_dx=float(np.max(np.abs(fx[1:] - (2 + (X[1:] - x0) * sc)))),
                         max_dy=float(np.max(np.abs(fy[1:] - (2 + (Y[1:] - y0) * sc)))),
                         bbox=[int(x0), int(x1), int(y0), int(y1)], scale=sc)

# runs of equal heading over steps n=2..N  ==  1, prime gaps, tail
h = head[2:]
brk = np.nonzero(h[1:] != h[:-1])[0] + 1
runs = np.diff(np.concatenate([[0], brk, [h.size]]))
gaps = np.diff(primes)
expect = np.concatenate([[1], gaps, [N - primes[-1]]]) if N > primes[-1] else np.concatenate([[1], gaps])
out["runs_equal_gaps"] = bool(runs.size == expect.size and np.all(runs == expect))
out["n_runs"] = int(runs.size)
vals, cnts = np.unique(gaps, return_counts=True)
order = np.argsort(-cnts)[:8]
out["gap_hist_top"] = [[int(vals[i]), int(cnts[i])] for i in order]
out["gap_max"] = int(gaps.max())

def path_stats(g):
    """g: run lengths (first run heading 0). Returns turning points, end, stats."""
    k = np.arange(g.size)
    dxs = np.array([1, 0, -1, 0])[k % 4] * g
    dys = np.array([0, -1, 0, 1])[k % 4] * g
    tx = np.concatenate([[0], np.cumsum(dxs)]); ty = np.concatenate([[0], np.cumsum(dys)])
    return tx, ty

def full_path(g):
    hd = np.repeat(np.arange(g.size) % 4, g)
    px = np.concatenate([[0], np.cumsum(np.array([1, 0, -1, 0])[hd])])
    py = np.concatenate([[0], np.cumsum(np.array([0, -1, 0, 1])[hd])])
    return px, py

def stats(g):
    tx, ty = path_stats(g)
    px, py = full_path(g)
    r = np.hypot(px, py)
    code = (px.astype(np.int64) + 10 ** 6) * 4 * 10 ** 6 + (py + 10 ** 6)
    distinct = np.unique(code).size
    tcode = (tx + 10 ** 6) * 4 * 10 ** 6 + (ty + 10 ** 6)
    tdist = np.unique(tcode).size
    rect = int(np.sum((g[:-3] == g[2:-1]) & (g[1:-2] == g[3:])))
    cls = np.array([g[c::4].sum() for c in range(4)])
    return dict(end_r=float(math.hypot(px[-1], py[-1])), max_r=float(r.max()),
                width=int(px.max() - px.min()), height=int(py.max() - py.min()),
                distinct_points=int(distinct), revisits=int(px.size - distinct),
                turning_point_repeats=int(tx.size - tdist), rect_loops=rect,
                class_spread=int(cls.max() - cls.min()))

g_real = runs.copy()
real = stats(g_real)
out["real"] = real
tx, ty = path_stats(g_real)
out["end"] = [int(X[N]), int(Y[N])]
cls = [int(g_real[c::4].sum()) for c in range(4)]
out["numbers_per_heading"] = cls   # heading c = pi(n-1) mod 4 ; plus n=1 at origin
out["primes_per_arrival_heading"] = [int(np.sum(head[primes] == c)) for c in range(4)]

# gap autocovariance of normalised gaps (gap / log p) and the turtle diffusion constant
gn = gaps / np.log(primes[:-1].astype(float))
gn = gn - gn.mean()
c = [float(np.mean(gn[: gn.size - l] * gn[l:])) for l in range(0, 9)]
out["autocov_norm_gaps"] = c
out["autocorr_norm_gaps"] = [x / c[0] for x in c]
D = c[0] + 2 * sum(c[l] * math.cos(math.pi * l / 2) for l in range(1, 9))
out["diffusion_ratio_D_over_c0"] = D / c[0]

# block-shuffled models
rng = np.random.default_rng(7)
NS = 400
rec = {k: [] for k in real}
B = 500
for s in range(NS):
    g = g_real.copy()
    core = g[1:-1]
    for a in range(0, core.size, B):
        blk = core[a:a + B]
        core[a:a + B] = blk[rng.permutation(blk.size)]
    g[1:-1] = core
    st = stats(g)
    for k in rec:
        rec[k].append(st[k])
summ = {}
for k, v in rec.items():
    v = np.array(v, float)
    summ[k] = dict(real=real[k], model_mean=float(v.mean()), model_sd=float(v.std()),
                   q05=float(np.quantile(v, .05)), q95=float(np.quantile(v, .95)),
                   frac_model_le_real=float(np.mean(v <= real[k])))
out["shuffle_model"] = summ
json.dump(out, open("turtle_out.json", "w"), indent=1)
for k, v in out.items():
    print(k, ":", v)
