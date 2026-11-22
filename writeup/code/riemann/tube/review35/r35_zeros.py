import mpmath as mp, math
mp.mp.dps = 20
chi = [0, 1, -1, -1, 1]; q = 5
def Z(t):
    s = mp.mpc(0.5, t)
    th = mp.im(mp.loggamma(s / 2)) + (t / 2) * mp.log(q / mp.pi)
    return mp.re(mp.exp(1j * th) * mp.dirichlet(s, chi))
T = 120; step = 0.04
ts = [i * step for i in range(1, int(T / step) + 1)]
vals = [Z(t) for t in ts]
zs = []
for i in range(len(ts) - 1):
    if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
        zs.append(mp.findroot(Z, (ts[i], ts[i + 1]), solver='anderson'))
zs = [float(z) for z in zs]
print(len(zs), 'zeros up to', T, 'first', zs[:4])
# Riemann-von Mangoldt count for even chi: N(T) ~ (T/pi) ln(qT/(2 pi e)) + O(...) (counting both signs: 2N+ ~ (T/pi)ln(qT/2pi e))
print('expected count of positive ordinates ~', (T / (2 * math.pi)) * math.log(q * T / (2 * math.pi * math.e)) + 0)
S = sum(2 / (0.25 + g * g) for g in zs)
# tail: density of positive ordinates (1/2pi) ln(q t/2pi); tail sum of 2/t^2 * density
tail = mp.quad(lambda t: 2 / t**2 * mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi), [T, mp.inf])
print('sum 1/|rho|^2 over found zeros', S, ' + tail', float(tail), '=', S + float(tail), ' vs C 0.1565570')
d1 = 0.8276794767155048879910469896773736467539; B0 = -0.02202465783872691842463476000992357551451
from sympy import primerange
def G(y): return math.fsum([math.log(p) * chi[p % 5] ** k * (y / p ** k - 1) for p in primerange(2, int(y) + 1) for k in range(1, 64) if p ** k <= y] + [y * d1])
def Te(y):
    u = 1 / y; return -0.5 * math.log1p(-u * u) - math.atanh(u) / u + 1
for y in [1.7, 2.5, 4, 7, 12, 30]:
    zsum = math.sqrt(y) * sum(2 * math.cos(g * math.log(y)) / (0.25 + g * g) for g in zs)
    lhs = G(y) - math.log(y) - 1 - B0 + Te(y)
    print(y, 'G - ln y - 1 - B0 + T_e =', round(lhs, 5), ' zero sum (truncated) =', round(zsum, 5), '  with -ln y instead:', round(G(y) + math.log(y) - 1 - B0 + Te(y), 5))
