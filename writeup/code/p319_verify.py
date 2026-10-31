# Puzzle 319: Morawetz's approximations to pi(N) are approximants of the asymptotic series of li(N).
import mpmath as mp, sympy as sp
mp.mp.dps = 40
u = sp.symbols('u'); A = 1 / u; B = u + 2*u**2 + 3*u**3 + 4*u**4 - 5*u**5
F = {'a': A * (sp.Rational(1, 2) - sp.sqrt(sp.Rational(1, 4) - u)),         # each divided by N/A
     'b': (A - 1) / (A - 2),
     'c': (A - 1 - sp.log(10) / A) / (A - 2 - sp.log(10) / A),
     'd': 1 + 1 / (A - 2 - B + sp.Rational(246, 10) / A**5)}             # our reading of the garbled (d)
co = {k: [sp.nsimplify(sp.series(f, u, 0, 6).removeO().coeff(u, i)) for i in range(6)] for k, f in F.items()}
assert co['a'] == [1, 1, 2, 5, 14, 42] and co['b'] == [1, 1, 2, 4, 8, 16] and co['d'][:5] == [1, 1, 2, 5, 14]
assert co['c'][:3] == [1, 1, 2] and abs(float(co['c'][3]) - (4 + float(sp.log(10)))) < 1e-12   # 6.3026 vs 3! = 6
# accuracy against exact pi(10^k) (OEIS A006880)
pis = {10: 455052511, 12: 37607912018, 14: 3204941750802, 16: 279238341033925, 18: 24739954287740860, 20: 2220819602560918840, 22: 201467286689315906290}
for k, pk in pis.items():
    N = mp.mpf(10)**k; a = mp.log(N); b = 1/a + 2/a**2 + 3/a**3 + 4/a**4 - 5/a**5
    c = (N/a) * (a - 1 - mp.log(10)/a) / (a - 2 - mp.log(10)/a)
    rel = float((c - pk) / pk)
    if k == 14: assert abs(rel - 1.06e-7) < 0.01e-7                   # Terr's figure on the page
    if k >= 16: assert 0.9e-6 < rel < 1.3e-6                           # (c) saturates near 1e-6
    assert abs(float((mp.li(N) - pk) / pk)) < 1e-5
print("all claims verified")
