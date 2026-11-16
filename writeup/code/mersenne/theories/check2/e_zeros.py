# Item E (zeros, V): L(s, (./P)) for primes P. Three routes for V.
import math, cypari2
import mpmath as mp
pari = cypari2.Pari()
pari.allocatemem(400 * 10**6)
mp.mp.dps = 30
PS = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
TZ = 400   # height for the zero-sum route

def chi_table(P):
    qr = {(a * a) % P for a in range(1, P)}
    return [0] + [1 if a in qr else -1 for a in range(1, P)]

print('P   D   a  first zero   V(Stieltjes/mpmath)  V(PARI lfun deriv)  V(zero sum to T=%d + tail)  #zeros' % TZ)
for P in PS:
    D = P if P % 4 == 1 else -P
    a = 0 if D > 0 else 1
    chi = chi_table(P)
    # sanity: Legendre symbol == Kronecker (D/n)
    for n in range(1, 200):
        if n % P:
            assert int(pari.kronecker(D, n)) == chi[n % P], (P, n)
    L = pari.lfuncreate(D)
    z = pari.lfunzeros(L, [0, 10])
    first = float(z[0])
    # route 1: L'/L(1) via Hurwitz expansion: L(s) = P^-s sum chi(a) zeta(s, a/P),
    # zeta(s,x) = 1/(s-1) - psi(x) - gamma_1(x)(s-1) + ...
    num = mp.fsum(chi[k] * mp.stieltjes(1, mp.mpf(k) / P) for k in range(1, P))
    den = mp.fsum(chi[k] * mp.digamma(mp.mpf(k) / P) for k in range(1, P))
    LL1 = -mp.log(P) + num / den
    V1 = 2 * (LL1 + mp.log(mp.mpf(P) / mp.pi) / 2 + mp.digamma(mp.mpf(1 + a) / 2) / 2)
    # route 2: PARI derivative
    l0 = pari.lfun(L, 1); l1 = pari.lfun(L, 1, 1)
    LL2 = float(l1 / l0)
    V2 = 2 * (LL2 + 0.5 * math.log(P / math.pi) + 0.5 * float(mp.digamma((1 + a) / 2)))
    # route 3: sum over zeros 1/|rho|^2 = 2 * sum_{gamma>0} 1/(1/4+gamma^2) + tail
    zs = [float(t) for t in pari.lfunzeros(L, [0, TZ])]
    s3 = 2 * sum(1 / (0.25 + g * g) for g in zs)
    # tail: zero density (1/pi) ln(P t / 2pi) dt for |gamma| (both signs together)
    tail = float(mp.quad(lambda t: mp.log(P * t / (2 * mp.pi)) / mp.pi / (0.25 + t * t), [TZ, mp.inf]))
    # correction for the count error at T: (N(T) - Nsmooth(T)) / (1/4 + T^2)
    Nsm = TZ / math.pi * math.log(P * TZ / (2 * math.pi * math.e)) + (0 if a == 0 else 0)
    V3 = s3 + tail
    print('%2d %4d %d  %.4f      %.5f              %.5f              %.5f                %d (N_smooth %.1f)' % (
        P, D, a, first, float(V1), V2, V3, 2 * len(zs), Nsm))
