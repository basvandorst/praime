"""Compare every C timing output speed_first_K_mode_filter.txt with the Python table ms_first.txt
(smallest k <= 1000; truncated to K).  Python twin of ms_speed.c = ms_main.py (numpy) and ms_check.py (pure Python)."""
import glob
ref = {}
for line in open("chk_first.txt"):
    p, k = map(int, line.split()); ref[p] = k
for fn in sorted(glob.glob("speed_first_*.txt")):
    K = int(fn.split("_")[2])
    got = {}
    for line in open(fn):
        p, k = map(int, line.split()); got[p] = k
    want = {p: (k if 0 < k <= K else 0) for p, k in ref.items()}
    bad = [p for p in want if want[p] != got.get(p)]
    print(f"{fn}: {'OK' if not bad else 'WRONG for %d primes, first %s' % (len(bad), bad[:3])}")
