"""Compare the re-polished optima (hp_d*.json, arb + own K0/K1) with the stored ones (certs/polish_A_d*.json,
mpmath.besselk): lambda, mu, touch points, coefficients b (relative differences)."""
import json, os
import common
from mpmath import mp, mpf
mp.dps = 130
print(" deg  |dlam|/lam   max|dmu|      max|du|      max|db_j|/|b_j|  stored-digits(b)  max|z|(new)")
for d in range(2, 26, 2):
    H = json.load(open("hp_d%d.json" % d)); S = json.load(open(os.path.join(common.CERTS, "polish_A_d%d.json" % d)))
    dl = abs(mpf(H["lam"]) - mpf(S["lam"])) / mpf(H["lam"])
    dm = max(abs(mpf(a) - mpf(b)) for a, b in zip(H["mu"], S["mu"]))
    du = max(abs(mpf(a) - mpf(b)) for a, b in zip(H["u"], S["touch_points"]))
    db = max(abs(mpf(a) - mpf(b)) / abs(mpf(a)) for a, b in zip(H["b"], S["b"]))
    print(" %3d  %.2e    %.2e    %.2e    %.2e          %d              %s" % (2 * d, dl, dm, du, db, len(S["b"][0]) - 2, max(H["z"], key=lambda x: abs(float(x)))))
