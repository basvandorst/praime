import numpy as np
from math import log

C2 = 0.6601618158468696
gaps = [2, 6, 12, 18, 30, 36, 72, 150, 168, 210, 282, 372, 498, 630, 924, 930, 1008, 1452, 1512, 1530, 1722, 1902, 2190, 2256, 2832, 2868, 3012, 3102, 3180, 3480]
starts = [3, 5, 17, 41, 71, 311, 347, 659, 2381, 5879, 13397, 18539, 24419, 62297, 187907, 687521, 688451, 850349, 2868959, 4869911, 9923987, 14656517, 17382479, 30752231, 32822369, 96894041, 136283429, 234966929, 248641037, 255949949]
g = np.array(gaps, dtype=float)
p1 = np.array(starts, dtype=float)
p2 = p1 + g
ok = p1 > 1000
lp1 = np.log(p1[ok])
lp2 = np.log(p2[ok])
rod = g[ok] / lp1 ** 3
a = lp2 ** 2 / (2 * C2)
kour = a * (np.log(p2[ok] / a) - 1.0)
print("record gaps between twin primes (OEIS A113274, A113275), lesser member of the first pair above 1000:")
print("  gap    p1            gap/ln^3 p1   Kourbatov a(ln(p2/a) - 1)   gap/formula")
for gg, pp, r, k in zip(g[ok], p1[ok], rod, g[ok] / kour):
    print(f"  {int(gg):5d}  {int(pp):12d}   {r:.3f}         {gg / k:10.0f}   {k:.3f}")
print(f"Rodriguez's ratio gap/ln^3 p1: min {rod.min():.3f}, max {rod.max():.3f}, mean {rod.mean():.3f} over {ok.sum()} records; it is not a constant but (1/(2C2)) (1 - (2 ln ln p + ln(2C2) + b)/ln p + ...) with limit 1/(2C2) = {1/(2*C2):.3f}")
print(f"Kourbatov's formula with b = 1: gap/formula mean {(g[ok]/kour).mean():.3f}, sd {(g[ok]/kour).std():.3f}")
print(f"ratio implied by Kourbatov's formula at p = 1e9, 1e12, 1e15, 1e20, 1e30: " + ", ".join(f"{(1/(2*C2))*(1 - (2*np.log(log(x)) - log(2*C2) + 1)/log(x)):.3f}" for x in (1e9, 1e12, 1e15, 1e20, 1e30)))
