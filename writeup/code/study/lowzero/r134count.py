"""Round 134, the registered report of the completeness repair: PARI lfunzeros at its default subdivision for every
discriminant of test B, counted against the saved (repaired) sets. Writes r134_repair.json; can stop and continue.
usage: r134count.py [SECONDS]"""
import sys, os, json, time
import numpy as np, cypari2
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros')
Z = np.load(os.path.join(ZD, 'r134_B.npz')); OUT = os.path.join(H, 'r134_repair.json')
log = json.load(open(OUT)) if os.path.exists(OUT) else {}
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6); budget = float(sys.argv[1]) if len(sys.argv) > 1 else 1e9; t0 = time.time()
for k in Z.files:
    if k in log: continue
    if time.time() - t0 > budget: break
    n0 = sum(1 for t in pari.lfunzeros(pari.lfuncreate(int(k[1:])), 600) if float(t) > 0)
    log[k] = [n0, int(len(Z[k]))]; json.dump(log, open(OUT, 'w'))
done = len(log); fixed = {k: v for k, v in log.items() if v[0] != v[1]}
print('%d of %d sets recounted; completed by the repair: %d sets, %d zeros added %s' % (done, len(Z.files), len(fixed), sum(v[1] - v[0] for v in fixed.values()), sorted(fixed)))
