import math, numpy as np
from sympy import primerange
X = 2*10**6
d1 = 0.8276794767155048879910469896773736467539; B0 = -0.02202465783872691842463476000992357551451
chi = [0, 1, -1, -1, 1]
ev = []
for p in primerange(2, X+1):
    v, k = p, 1
    while v <= X:
        ev.append((v, math.log(p) * chi[p % 5] ** k)); v *= p; k += 1
ev.sort()
n = np.array([e[0] for e in ev], float); a = np.array([e[1] for e in ev])
A = np.cumsum(a / n); B = np.cumsum(a)
al = A + d1; be = -B - 1 - B0
xr = np.append(n[1:], X)
h = lambda t: al * t - be + np.log(t) - 2
hl, hr = h(n), h(xr)
m = n >= 100
print('sign changes inside intervals (n>=100):', int(np.sum((np.sign(hl) != np.sign(hr)) & m)))
print('min |h| at left ends:', np.min(np.abs(hl[m])), ' typical |h| median', np.median(np.abs(hl[m])))
# finite difference check of g' sign vs h at random interior points
rng = np.random.default_rng(1); idx = rng.choice(np.nonzero(m)[0], 5)
for i in idx:
    t = (n[i] + xr[i]) / 2; g = lambda y: (al[i]*y + be[i] - math.log(y)) / math.sqrt(y)
    print(n[i], (g(t+1e-4)-g(t-1e-4))/2e-4, h(t)[i] / (2*t**1.5))
# frozen-ln version of the repo
cc = np.log(n) + 1 + B0; beta = -(B + cc); xs = beta / al
ok = (xs > n) & (xs < xr) & (xs >= 100)
print('repo frozen-ln interior candidates:', int(ok.sum()))
if ok.any():
    t = xs[ok]; rv = np.abs(al[ok]*t - B[ok] - (np.log(t)+1+B0))/np.sqrt(t); print(t[:10], rv[:10])
