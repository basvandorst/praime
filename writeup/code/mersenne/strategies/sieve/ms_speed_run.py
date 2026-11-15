"""Run ms_speed (C) for a set of (K, mode, filter), 3 repeats each, report the fastest; then ms_speed_check.py validates outputs."""
import subprocess, re, os
cfg = [(64, m, f) for m in (0, 1, 3, 4) for f in (0, 1, 2)] + [(94, 1, f) for f in (0, 1, 2)] + \
      [(256, m, f) for m in (2, 3, 4) for f in (0, 1, 2)] + [(1000, m, 2) for m in (3, 4)]
print("load:", os.getloadavg())
for K, m, f in cfg:
    best = None
    for _ in range(3):
        out = subprocess.run(["./ms_speed", str(K), str(m), str(f)], capture_output=True, text=True).stdout
        t = float(re.search(r"([\d.]+) s,", out).group(1)); mp = int(re.search(r"(\d+) modpows", out).group(1))
        best = t if best is None else min(best, t)
    print(f"K={K:5d} mode={m} filter={f}: modpows {mp:9d}  best of 3: {best:.3f} s  ({1e9*best/mp:.0f} ns/modpow)")
print("load:", os.getloadavg())
