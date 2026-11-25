"""Our level formula against BBLM's at the heights of their figures (E = 2.504e15: L = 33.62; E = 1.307e22: L = 49.09)
and at L = 20, 100: Var(s) deviations from GUE and the gap ratio."""
import sys, os
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/lib'); sys.path.insert(0, H)
import level
lam, Q, V0, R0 = 1.573151, 2.315846, level.V0, level.R0
print('   L      dVar ours   BBLM Lambda   BBLM Lambda+Q  |  d<r~> ours   (N_eff, alpha)')
for L in (20.0, 33.62, 49.09, 100.0):
    dv, dr = level.level(L)
    print('  %6.2f  %+.6f  %+.6f  %+.6f  |  %+.6f   (%.3f, %.4f)' % (L, dv, -(2 * lam / L ** 2) * (1 + V0), -(2 * lam / L ** 2 - 6 * Q / L ** 3) * (1 + V0), dr, L / (12 * lam) ** 0.5, 1 + Q / (lam * L)))
