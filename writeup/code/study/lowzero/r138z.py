"""Round 138: zeros (0 < gamma <= TOP) of the L-functions of the fundamental discriminants of test Fn or Bp
(triple/predictions_kappa3.json), PARI lfunzeros, then the completeness repair of fixpairs.py (count drift against
N(t) above 1.2: the 30 largest unfolded gaps are searched again with 256 subdivisions). One file per discriminant, so a
run can stop and continue. usage: r138z.py Fn|Bp PART NPARTS [SECONDS]; r138z.py Fn|Bp pack"""
import sys, os, json, math, time
import numpy as np, cypari2
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros')
T = sys.argv[1]; P = json.load(open(os.path.join(H, 'triple', 'predictions_kappa3.json')))[T]
Ds = P['D']; TOP = {'Fn': 600, 'Bp': 1000}[T]; OUT = os.path.join(ZD, 'r138_' + T); os.makedirs(OUT, exist_ok=True)
def drift(g, q):
    g = np.sort(g[g > 20]); u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); r = np.arange(len(g)) - (u - u[0])
    k = max(20, len(g) // 20); return float(r[-k:].mean() - r[:k].mean())
if sys.argv[2] == 'pack':
    Z = {'D%d' % D: np.load(os.path.join(OUT, 'D%d.npy' % D)) for D in Ds}; log = json.load(open(os.path.join(OUT, 'log.json'))) if os.path.exists(os.path.join(OUT, 'log.json')) else {}
    np.savez_compressed(os.path.join(ZD, 'r138_%s.npz' % T), **Z)
    dr = {D: drift(Z['D%d' % D], abs(D)) for D in Ds}
    print(len(Z), 'sets packed;', sum(len(v) for v in Z.values()), 'zeros; |drift| > 1.2 in', sum(abs(v) > 1.2 for v in dr.values()), 'sets; max |drift| %.2f' % max(abs(v) for v in dr.values()))
    sys.exit()
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
part, nparts = int(sys.argv[2]), int(sys.argv[3]); budget = float(sys.argv[4]) if len(sys.argv) > 4 else 1e9; t0 = time.time()
def fix(g, D):
    q = abs(D); added = 0
    for _ in range(3):
        if abs(drift(g, q)) <= 1.2: break
        gs = np.sort(g); u = gs / (2 * math.pi) * np.log(q * gs / (2 * math.pi * math.e)); z = []
        for i in np.argsort(-np.diff(u))[:30]:
            z += [float(t) for t in pari('lfunzeros(lfuncreate(%d),[%f,%f],256)' % (D, gs[i] + 1e-6, gs[i + 1] - 1e-6))]
            if z: break
        if not z: break
        g = np.sort(np.concatenate([g, z])); added += len(z)
    return g, added
done = 0
for D in Ds[part::nparts]:
    f = os.path.join(OUT, 'D%d.npy' % D)
    if os.path.exists(f): continue
    if time.time() - t0 > budget: break
    g = np.array(sorted(float(t) for t in pari.lfunzeros(pari.lfuncreate(D), TOP) if float(t) > 0)); n0 = len(g)
    g, added = fix(g, D); np.save(f, g); done += 1
    print(D, n0, added, '%.2f' % drift(g, abs(D)), '%.0fs' % (time.time() - t0), flush=True)
print('part', part, 'computed', done, 'left', sum(not os.path.exists(os.path.join(OUT, 'D%d.npy' % D)) for D in Ds[part::nparts]))
