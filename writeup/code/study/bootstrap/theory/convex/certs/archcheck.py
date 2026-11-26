"""Independent check of every polished optimum: cost by the prime sum (closed form, core.Basis, recomputed here
from scratch at DPS (default 80) digits, prime powers n <= 400) against the archimedean formula
   cost = B(i/2)/2 + sum_j b_j a_j,   a_j = (1/2pi) int Xi^2 t^{2j} (Re psi(1/4+it/2) - log pi) dt,
and g(0) = sum_j b_j m_j (Xi side) against g(0) = sum_j b_j g_j(0) (closed form).
Moments from ximom.py (file XIMOM, default xi_mom48.pkl: 48 Gauss-Legendre nodes per unit panel; xi_mom64.pkl: 64 nodes, 110 digits).
Usage: python archcheck.py"""
import json, os, glob, pickle
import mpmath as mp
from core import Basis, von_mangoldt

mp.mp.dps = int(os.environ.get("DPS", 80))
XM = pickle.load(open(os.environ.get("XIMOM", "xi_mom48.pkl"), "rb"))
aj = [mp.mpf(x) for x in XM["a"]]
mj = [mp.mpf(x) for x in XM["m"]]
print("Xi moments: J=%d dps=%d T=%d nodes/panel=%s" % (XM["J"], XM["dps"], XM["T"], XM.get("nodes", 24)))
files = sorted(glob.glob("polish_*_d*.json"), key=lambda s: (s.split("_")[1], int(s.split("_d")[1][:-5])))
Jmax = max(len(json.load(open(f))["b"]) - 1 for f in files)
Jmax = min(Jmax, XM["J"])
bas = Basis(Jmax, 0)
PPS = [m for m in range(2, 401) if von_mangoldt(m)]
vals = {m: [x[0] for x in bas.eval(mp.log(m))] for m in PPS}
v0 = [x[0] for x in bas.eval(0)]
W = {m: 2 * von_mangoldt(m) / mp.sqrt(m) for m in PPS}
for f in files:
    r = json.load(open(f))
    b = [mp.mpf(x) for x in r["b"]]
    n = len(b)
    if n - 1 > Jmax:
        print("%-20s skipped (needs moments up to j=%d)" % (f, n - 1))
        continue
    prime = mp.fsum(W[m] * mp.fsum(b[j] * vals[m][j] for j in range(n)) for m in PPS)
    g0 = mp.fsum(b[j] * v0[j] for j in range(n))
    arch = mp.fsum(b[j] * ((-mp.mpf(1) / 4) ** j / 2 + aj[j]) for j in range(n))
    g0xi = mp.fsum(b[j] * mj[j] for j in range(n))
    print("%-18s prime: %s  arch: %s  rel.diff %s  | g0 closed-Xi %s" % (
        f, mp.nstr(prime / g0, 16), mp.nstr(arch / g0xi, 16), mp.nstr((prime / g0 - arch / g0xi) / (prime / g0), 3),
        mp.nstr(g0 - g0xi, 3)), flush=True)
