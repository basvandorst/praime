# Problem 88 Q3: K with 8aK^2+1 and 8bK^2+1 both squares (a = 222121, b = 360721), K > 0.
# K = k_n(a) = k_m(b), k_n = (eps^n - eps^-n)/(2 sqrt(8s)). Equality forces |n La - m Lb - (1/2)ln(a/b)| < 4 eps_min^(-2).
from sympy.solvers.diophantine.diophantine import diop_DN
from mpmath import mp, mpf, log, sqrt, nint, fabs
a, b = 222121, 360721
fa = min(diop_DN(8 * a, 1), key=lambda t: t[0]); fb = min(diop_DN(8 * b, 1), key=lambda t: t[0])
print("fundamental (h, k):", fa, fb)
assert fa == (22632419659503, 16978181966)
assert fb == (86286108053753141137230600854902563, 50793737929709259543442603308701)
mp.dps = 80
La = log(fa[0] + fa[1] * sqrt(8 * a)); Lb = log(fb[0] + fb[1] * sqrt(8 * b)); c = log(mpf(a) / b) / 2
# common term with n, m >= 1 needs |n La - m Lb - c| < 4 * min(eps)^(-2) ~ 1e-27
best = mpf(1); arg = None
N = 10**6
for n in range(1, N + 1):
    m = nint((n * La - c) / Lb)
    if m < 1: continue
    d = fabs(n * La - m * Lb - c)
    if d < best: best, arg = d, (n, int(m))
print(f"min over 1 <= n <= {N} of |n La - m Lb - c| = {mp.nstr(best, 5)} at (n, m) = {arg}")
print("required for a solution: < 4*eps_a^-2 =", mp.nstr(4 / (fa[0] + fa[1] * sqrt(8 * a))**2, 5))
