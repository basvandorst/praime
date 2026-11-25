import math, cmath, numpy as np, sys
sys.argv = ['x']
src = open('r3_g13.py').read().split('# ---------------------------------------------------------------- chi mod 7')[0]
exec(src)
L7 = P.lfuncreate(P('[znstar(7,1),[1]]')); L7b = P.lfuncreate(P('[znstar(7,1),[5]]'))
chi = np.zeros(7, complex)
for n in range(1, 7): chi[n] = cmath.exp(2j * PI * float(P(f'chareval(znstar(7,1),[1],{n})')))
z7 = np.array([float(x) for x in P.lfunzeros(L7, P('[0,60]'))]); z7b = np.array([float(x) for x in P.lfunzeros(L7b, P('[0,60]'))])
print('lowest zeros', z7[:3], z7b[:3])
for N in (40000, 400000):
  pr = primes(N); ns, Lam = [], []
  for p in pr:
    p = int(p); q = p; k = 1
    while q <= N: ns.append(q); Lam.append(chi[p % 7] ** k * math.log(p)); q *= p; k += 1
  ns = np.array(ns, float); Lam = np.array(Lam)
  for h in (0.05, 0.02):
   for cut in (6.5, 9.0):
    T1, T2, D, r = 3.0, 15.0, 1.5, 1.0
    a, b = T1 - cut * D, T2 + cut * D
    t = a + h * np.arange(int((b - a) / h) + 1)
    f = Hbox(t, T1, T2, D) * e2it(t, 7.0, [1.0])
    g = np.concatenate([z7, -z7b])
    Z = np.sum(Hbox(g, T1, T2, D) * e2it(g, 7.0, [1.0]))
    I = horner(np.log(ns), t[0], h, f); Ps = -np.sum(Lam * ns ** -0.5 * I) / (2 * PI)
    thp, _ = thetap(t, 7.0, [1.0])
    Im_ = horner(np.log(ns), t[0], h, f, sign=+1)
    E0 = h * np.sum(f * thp) / PI - np.sum(np.conj(Lam) * ns ** -0.5 * Im_) / (2 * PI)
    print(f'N={N} h={h} cut={cut}: |Z-P-E0| = {abs(Z-Ps-E0):.2e}, H(a)={Hbox(a,T1,T2,D):.1e}')
