#!/usr/bin/env python3
"""rh_density_check.py -- numeric checks for rh_density.md (zero-density route to Theorem 1).

Own code; imports only numpy, mpmath and the standard library.
Run:  (ulimit -v 1000000; timeout 600 nice python -I -B rh_density_check.py)

Sections
  A  exponent of the zero-detection method as a function of kappa; explicit theta values (exact fractions)
  B  constants of the classical lemmas (digamma, Binet, Gamma decay, zero-count constants, tails)
  C  functional equation of Hecke L-functions of Q(i) with Groessencharacters: the Gamma factor
     Gamma(s+|k|/2) and A = sqrt(4 Nf)/(2 pi) (root number |W| = 1, independent of s and of the split);
     negative controls; sign changes on the critical line vs the zero-count bound of Lemma ZC
  D  zero-detection identities and inequalities on real and synthetic Dirichlet polynomials;
     Lemma 5 of rh_proof.md on the coefficient vectors that the density proof uses
  E  ideal arithmetic of Z[i] used in the bookkeeping
  F  sanity only: Gaussian primes in classes mod 2+i and mod 5 in sectors, N z <= 1e7; sizes of Psi(xi)
  G  illustration (heuristic, placeholder constants): VK region vs classical region in the zero sum
  H  the inequalities repaired after the referee report: consequence (a) of Lemma PF (neighbouring windows);
     (T1)-(T3) of section 4.3 (divisor sums, r0 <= d, the dyadic tail bound for all Y >= 100); the bound (Gamma)
     of Lemma DI on -1/2 <= Re w <= 2 through the shift Gamma(w) = Gamma(w+1)/w, and (C3)(iii); the block sum of 5
"""
import math, cmath, random, time
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp

T_START = time.time()
NCHK = 0
NFAIL = 0


def check(name, ok, info=""):
    global NCHK, NFAIL
    NCHK += 1
    if not ok:
        NFAIL += 1
    print(("PASS " if ok else "FAIL ") + name + ((" | " + info) if info else ""))


def hdr(s):
    print("\n=== " + s + " ===")


# ---------------------------------------------------------------- A. exponents
hdr("A. exponent of the method")


def detection_exponent(a, b, sig, kap):
    """log_M of the zero-detection bound (logs dropped) for X = M^a, Y = M^b, vectorised over sigma.
    type I : levels U = M^u, u in [a, b]; for u >= 1: (M+U)U^(1-2s) -> u(2-2s); for u < 1 the best of
             powering by g = 1, 2, ...: max(1, g u) + g u (1-2s)  (g = 1 is no powering);
    type II: (M+X) Y^(1-2s) M^(4 kappa) -> max(1,a) + 4 kappa + b(1-2s)."""
    s = np.asarray(sig)[:, None, None]
    u = np.linspace(a, b, 61)[None, :, None]
    g = np.arange(1, 12)[None, None, :]
    pw = np.maximum(1.0, g * u) + g * u * (1 - 2 * s)
    small = pw.min(axis=2)
    big = (u * (2 - 2 * s))[:, :, 0]
    t1 = np.where(u[:, :, 0] >= 1, big, small).max(axis=1)
    t2 = max(1.0, a) + 4 * kap + b * (1 - 2 * np.asarray(sig))
    return np.maximum(t1, t2)


def method_A(a, b, kap, sigs):
    e = np.minimum(1.0, detection_exponent(a, b, sigs, kap))   # trivial bound N << M T0 log: exponent 1
    return float((e / (1 - sigs)).max())


sigs = np.concatenate([np.linspace(0.5, 0.99, 99), 1 - np.logspace(-2, -5, 13)])
for kap in [0.02, 0.06, 0.125, 1 / 6, 0.2, 103 / 512, 0.24, 0.25]:
    A_choice = method_A(1.0, 1 + 4 * kap, kap, sigs)
    best = 1e9
    arg = None
    for a in np.linspace(0.5, 1.5, 11):
        for b in sorted(set(list(np.linspace(1.0, 2.5, 31)) + [1 + 4 * kap])):
            v = method_A(a, b, kap, sigs)
            if v < best:
                best, arg = v, (a, b)
    check("A(kappa=%.4f): choice X=M, Y=M^(1+4k) gives 2+8k" % kap, abs(A_choice - (2 + 8 * kap)) < 1e-9,
          "sup ratio %.9f vs %.9f" % (A_choice, 2 + 8 * kap))
    check("A(kappa=%.4f): grid min over (X,Y) not below 2+8k" % kap, best >= 2 + 8 * kap - 1e-9,
          "grid min %.6f at (a,b)=(%.2f,%.3f)" % (best, arg[0], arg[1]))

# both zero types equal at the choice
for kap in [1 / 6, 103 / 512]:
    for s in [0.7, 0.8, 0.95]:
        b = 1 + 4 * kap
        t1 = b * (2 - 2 * s)
        t2 = 1 + 4 * kap + b * (1 - 2 * s)
        check("types I and II equal at kappa=%.4f sigma=%.2f" % (kap, s), abs(t1 - t2) < 1e-12 and abs(t1 - (2 + 8 * kap) * (1 - s)) < 1e-12)

# exact fractions
def kappa_wu(vt):
    return Fr(1, 4) - (1 - 2 * vt) / 16


cases = [("Wu + Blomer-Brumley theta=7/64", kappa_wu(Fr(7, 64)), Fr(64, 231)),
         ("Wu + Kim-Shahidi theta=1/9", kappa_wu(Fr(1, 9)), Fr(18, 65)),
         ("Wu + Luo-Rudnick-Sarnak theta=3/10", kappa_wu(Fr(3, 10)), Fr(5, 19)),
         ("Wu + Ramanujan theta=0", kappa_wu(Fr(0)), Fr(2, 7)),
         ("Soehne kappa=1/6", Fr(1, 6), Fr(3, 10)),
         ("convexity kappa=1/4", Fr(1, 4), Fr(1, 4))]
for name, k, th in cases:
    A = 2 + 8 * k
    check("theta* = 1/(2+8k) for " + name, 1 / A == th, "kappa=%s A=%s theta*=%s=%.5f sigma1=1-1/A=%s" % (k, A, 1 / A, float(1 / A), 1 - 1 / A))
check("kappa(Wu, 7/64) = 103/512", kappa_wu(Fr(7, 64)) == Fr(103, 512))
check("64/231 > 1/4 and 3/10 > 64/231", Fr(64, 231) > Fr(1, 4) and Fr(3, 10) > Fr(64, 231))
# admissible parameters for one theta below 64/231
th, kap = 0.2770, 103 / 512 + 1e-6
eps = (1 / (2 + 8 * kap) - th) / 3
check("parameters (E): (theta+2eps)(2+8kappa) < 1 and theta+2eps < 1/2", (th + 2 * eps) * (2 + 8 * kap) < 1 and th + 2 * eps < 0.5,
      "theta=%.4f eps=%.2e delta'=%.2e" % (th, eps, 1 - (th + 2 * eps) * (2 + 8 * kap)))
# sigma1 - 1/2 = 2k/(1+4k)
for k in [Fr(1, 6), Fr(103, 512)]:
    A = 2 + 8 * k
    check("sigma1 - 1/2 = 2k/(1+4k) for k=%s" % k, (1 - 1 / A) - Fr(1, 2) == 2 * k / (1 + 4 * k))

# ---------------------------------------------------------------- B. constants
hdr("B. constants of the classical lemmas")
mp.mp.dps = 30
# Binet function g(t) = 1/(e^t-1) - 1/t + 1/2 satisfies 0 <= g(t) <= t/12
with mp.workdps(60):
    tt = [mp.mpf(10) ** mp.mpf(e) for e in np.linspace(-6, 2, 400)]
    gv = [1 / mp.expm1(t) - 1 / t + mp.mpf(1) / 2 for t in tt]
    okb = all(0 <= g <= t / 12 for g, t in zip(gv, tt))
check("Binet g: 0 <= g(t) <= t/12 on (1e-6, 100) (60 digits, expm1)", okb)
# |psi(z) - log z| <= 4/3 for Re z >= 1/2
worst = 0
for x in [0.5, 0.6, 0.8, 1, 1.5, 2, 3, 5, 10, 30, 100]:
    for y in [0, 0.1, 0.5, 1, 2, 5, 10, 30, 100, 1000]:
        z = mp.mpc(x, y)
        worst = max(worst, abs(mp.digamma(z) - mp.log(z)))
check("|psi(z) - log z| <= 4/3 for Re z >= 1/2 (grid)", worst <= 4 / 3, "max %.4f" % worst)
# Lemma S(c): |psi(z)| <= log(|z|+3) + 8 for Re z >= -1/4, |z| >= 1/4
worst = -1e9
for x0 in np.linspace(-0.25, 3, 27):
    for y0 in np.concatenate([np.linspace(-3, 3, 61), np.array([-100, -30, -10, 10, 30, 100])]):
        z = mp.mpc(x0, y0)
        if abs(z) < 0.25:
            continue
        worst = max(worst, abs(mp.digamma(z)) - mp.log(abs(z) + 3))
check("|psi(z)| <= log(|z|+3) + 8 for Re z >= -1/4, |z| >= 1/4 (grid)", worst <= 8, "max |psi| - log(|z|+3) = %.3f" % worst)
# |mu(z)| <= 1/(12 Re z) (Stirling remainder)
worst = 0
for x in [0.5, 0.7, 1, 2, 5]:
    for y in [0, 0.5, 2, 10, 50]:
        z = mp.mpc(x, y)
        mu = mp.loggamma(z) - ((z - 0.5) * mp.log(z) - z + mp.log(2 * mp.pi) / 2)
        worst = max(worst, abs(mu) * 12 * x)
check("Stirling remainder |mu(z)| <= 1/(12 Re z)", worst <= 1 + 1e-12, "max 12 Re z |mu| = %.4f" % worst)
# |Gamma(a+iv)| <= 3 e^{-|v|} for 1/2 <= a <= 1
worst = 0
for a in np.linspace(0.5, 1, 11):
    for v in np.linspace(-60, 60, 1201):
        worst = max(worst, abs(mp.gamma(mp.mpc(a, v))) * mp.e ** abs(v))
check("|Gamma(a+iv)| <= 3 e^{-|v|}, 1/2<=a<=1", worst <= 3, "sup %.4f (proof gives e^1.086 = 2.96)" % worst)
# Re 1/(s1-rho) >= 2/5 for s1 = 2+iT, 0<=beta<=1, |gamma-T|<=1
m = min((2 - b) / ((2 - b) ** 2 + d ** 2) for b in np.linspace(0, 1, 101) for d in np.linspace(0, 1, 101))
check("Re 1/(s1-rho) >= 2/5 in the unit box", m >= 0.4 - 1e-15, "min %.6f" % m)


def zetaK_logder(s):
    # -zeta_K'/zeta_K(s) = -zeta'/zeta(s) - L'/L(s, chi_-4) for K = Q(i), s real > 1
    chi = [0, 1, 0, -1]
    z = -mp.diff(lambda u: mp.log(mp.zeta(u)), s)
    l = -mp.diff(lambda u: mp.log(mp.dirichlet(u, chi)), s)
    return z + l


c2 = zetaK_logder(2)
c54 = zetaK_logder(mp.mpf(5) / 4)
check("-zeta_K'/zeta_K(2) <= 1.14", c2 <= 1.14, "value %.6f" % c2)
check("-zeta_K'/zeta_K(5/4) <= 10", c54 <= 10, "value %.6f" % c54)
# Lemma ZC constant: 5/2 * (c2 - log pi + 4/3 + log 2-shift) <= 3.35 (see write-up)
cz = 2.5 * (float(c2) - math.log(math.pi) + 4 / 3)
check("Lemma ZC additive constant 5/2(c_K - log pi + 4/3) <= 3.35", cz <= 3.35, "value %.4f" % cz)
# integrals and sums
I1 = mp.quad(lambda v: (1 + abs(v)) * mp.e ** (-abs(v)), [-mp.inf, 0, mp.inf])
check("int (1+|v|) e^{-|v|} dv = 4", abs(I1 - 4) < 1e-20)
I2 = mp.quad(lambda T: 0.5 * mp.e ** (-T / 2) * (T + 3), [0, mp.inf])
check("int_0^inf (1/2) e^{-T/2}(T+3) dT = 5", abs(I2 - 5) < 1e-20)
sd = 2 + 2 / (math.e - 1)
check("sum_n e^{-d_n} = 2 + 2/(e-1) <= 3.2", sd <= 3.2, "%.4f" % sd)
cps = 2 + math.log(2) ** 2 / 2
check("partial-summation constant 2 + (log 2)^2/2 <= 2.25", cps <= 2.25, "%.4f" % cps)
check("e^{-1/Y} >= 0.99 for Y >= 100", math.exp(-1 / 100) >= 0.99)


def tail_bound(Y):
    # sum_{N a > Y1} d(a) N a^{-1/2} e^{-N a/Y} <= sum_j e^{-2^j Y1/Y} (2^j Y1)^{-1/2} 2^{j+1} Y1 (1+log(2^{j+1}Y1))^3
    Y1 = Y * math.log(Y) ** 2
    tot = 0.0
    for j in range(0, 200):
        u = 2 ** j * Y1
        term = math.exp(-u / Y) * u ** -0.5 * 2 * u * (1 + math.log(2 * u)) ** 3
        tot += term
        if term < 1e-300:
            break
    return tot


tb = max(tail_bound(Y) for Y in np.logspace(2, 12, 60))
check("mollified tail beyond Y1 = Y log^2 Y is <= 1e-4 for Y >= 100", tb <= 1e-4, "max %.2e" % tb)
okd = all((math.log2(Y * math.log(Y) ** 2 / X) + 2) <= 2 * math.log(Y * math.log(Y) ** 2)
          for Y in np.logspace(2, 12, 30) for X in [2, 10, math.sqrt(Y)])
check("number of dyadic levels <= 2 log Y1", okd)

# ---------------------------------------------------------------- C. functional equation, zero counts
hdr("C. Hecke L-functions of Q(i): Gamma factor, root number, zero counts")
mp.mp.dps = 25


def ideal_gens(Nmax):
    out = []
    R = math.isqrt(Nmax) + 1
    for a in range(1, R + 1):
        for b in range(0, R + 1):
            n = a * a + b * b
            if n <= Nmax:
                out.append((a, b, n))
    return out   # one generator (a>0, b>=0) per nonzero ideal


LOG5 = {1: 0, 2: 1, 4: 2, 3: 3}   # 2 generates (Z/5)^*


def chi_mod_2pi(a, b, j):
    """character chi_j of (Z[i]/(2+i))^* = (Z/5)^* via a+bi -> a-2b mod 5, chi_j(2)=i^j; chi_j(i)=i^(3j)."""
    r = (a - 2 * b) % 5
    if r == 0:
        return 0
    return 1j ** (j * LOG5[r])


def xi_val(a, b, k, j, mod):
    e = cmath.exp(1j * k * math.atan2(b, a))
    return e if mod == 1 else chi_mod_2pi(a, b, j) * e


def split_parts(s, c, k, j, mod, ids, gam_shift=None):
    Nf = 1 if mod == 1 else 5
    A = mp.sqrt(4 * Nf) / (2 * mp.pi)
    ak = abs(k) / mp.mpf(2) if gam_shift is None else gam_shift
    P = mp.mpc(0)
    Q = mp.mpc(0)
    for (a, b, n) in ids:
        v = xi_val(a, b, k, j, mod)
        if v == 0:
            continue
        y = n / A
        P += v * y ** (-s) * mp.gammainc(s + ak, c * y)
        Q += v.conjugate() * y ** (s - 1) * mp.gammainc(1 - s + ak, y / c)
    return P, Q


def root_number(s, k, j, mod, ids, c1=1.0, c2=0.6, gam_shift=None):
    P1, Q1 = split_parts(s, c1, k, j, mod, ids, gam_shift)
    P2, Q2 = split_parts(s, c2, k, j, mod, ids, gam_shift)
    return (P1 - P2) / (Q2 - Q1), P1, Q1


ids300 = ideal_gens(300)
# chi_j(i) = i^{-k_chi} with k_chi = j (mod 4): check
for j in range(4):
    ci = chi_mod_2pi(0, 1, j)
    check("chi_%d(i) = i^(-k_chi), k_chi = %d" % (j, j), abs(ci - 1j ** (-j)) < 1e-12)
Wtab = {}
for (mod, k, j) in [(1, 4, 0), (1, 8, 0), (1, -12, 0), (5, 1, 1), (5, 2, 2), (5, 3, 3), (5, 5, 1), (5, -1, 3), (5, 9, 1), (5, 6, 2)]:
    Ws = []
    for s in [mp.mpc(0.3, 1.0), mp.mpc(0.8, -4.0), mp.mpc(-0.5, 7.0)]:
        W, _, _ = root_number(s, k, j, mod, ids300)
        Ws.append(W)
    W2, _, _ = root_number(mp.mpc(0.3, 1.0), k, j, mod, ids300, c1=1.3, c2=0.8)
    spread = max(abs(w - Ws[0]) for w in Ws + [W2])
    Wtab[(mod, k, j)] = complex(Ws[0])
    check("FE with Gamma(s+|k|/2): modulus %s, k=%d: |W|=1, W independent of s and split" % ("1" if mod == 1 else "2+i", k),
          abs(abs(Ws[0]) - 1) < 1e-12 and spread < 1e-12, "W=%s, spread %.1e" % (mp.nstr(Ws[0], 8), float(spread)))
# negative controls: wrong Gamma factor; a pair (chi, k) with k != k_chi mod 4 (not a character on ideals)
Wbad = [root_number(s, 4, 0, 1, ids300, gam_shift=mp.mpf(4))[0] for s in [mp.mpc(0.3, 1.0), mp.mpc(0.8, -4.0)]]
check("negative control: Gamma(s+|k|) instead of Gamma(s+|k|/2) is detected", abs(abs(Wbad[0]) - 1) > 1e-3 or abs(Wbad[0] - Wbad[1]) > 1e-3,
      "W=%s, %s" % (mp.nstr(Wbad[0], 6), mp.nstr(Wbad[1], 6)))
Wbad2 = [root_number(s, 0, 2, 5, ids300)[0] for s in [mp.mpc(0.3, 1.0), mp.mpc(0.8, -4.0)]]
check("negative control: k=0 with chi_2 (k != k_chi mod 4) fails the FE", abs(abs(Wbad2[0]) - 1) > 1e-3 or abs(Wbad2[0] - Wbad2[1]) > 1e-3)
Wbad3 = [root_number(s, 4, 0, 5, ids300)[0] for s in [mp.mpc(0.3, 1.0), mp.mpc(0.8, -4.0)]]
check("negative control: imprimitive xi (chi_0 mod 2+i, k=4) fails the FE with Nf = 5 (primitive version has Nf = 1)",
      abs(abs(Wbad3[0]) - 1) > 1e-3 or abs(Wbad3[0] - Wbad3[1]) > 1e-3)

# sign changes on the critical line vs Lemma ZC bound (lower bound for zeros vs upper bound)
mp.mp.dps = 18
ids90 = ideal_gens(90)


def zc_bound(Nf, k, T):
    return 1.25 * math.log(Nf) + 2.5 * math.log(abs(k) + abs(T) + 2) + 3.35


for (mod, k, j) in [(1, 4, 0), (5, 1, 1)]:
    W = Wtab[(mod, k, j)]
    sw = cmath.sqrt(W)
    ts = [0.05 * i for i in range(0, 281)]
    vals = []
    imax = 0.0
    for t in ts:
        P, Q = split_parts(mp.mpc(0.5, t), 1.0, k, j, mod, ids90)
        v = complex((P + W * Q) / sw)
        vals.append(v.real)
        imax = max(imax, abs(v.imag) / (abs(v) + 1e-300))
    sc = [0.5 * (ts[i] + ts[i + 1]) for i in range(len(ts) - 1) if (vals[i] > 0) != (vals[i + 1] > 0)]
    worst = 0.0
    for T in np.linspace(0, 14, 57):
        cnt = sum(1 for g in sc if abs(g - T) <= 1)
        worst = max(worst, cnt / zc_bound(1 if mod == 1 else 5, k, T))
    check("critical line, modulus %s k=%d: W^{-1/2}Lambda real; sign changes per window <= ZC bound" % ("1" if mod == 1 else "2+i", k),
          imax < 1e-8 and worst <= 1, "%d sign changes in [0,14]: %s; max count/bound %.2f" % (len(sc), [round(g, 2) for g in sc], worst))

# ---------------------------------------------------------------- D. zero detection
hdr("D. zero-detection identities and inequalities")
mp.mp.dps = 20
chi4 = [0, 1, 0, -1]
rho1 = mp.findroot(lambda s: mp.dirichlet(s, chi4), mp.mpc(0.5, 6.02))
check("first zero of L(s, chi_-4) found on Re s = 1/2", abs(rho1.real - 0.5) < 1e-12 and abs(mp.dirichlet(rho1, chi4)) < 1e-15,
      "rho1 = %s" % mp.nstr(rho1, 15))


def mobius_int(n):
    r, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            r = -r
        p += 1
    return -r if n > 1 else r


Xm, Ym = 10, 30
mu_tab = [0] + [mobius_int(n) for n in range(1, 2001)]
cnum = [0] * 2001
for n in range(1, 2001):
    cnum[n] = sum(mu_tab[d] for d in range(1, Xm + 1) if n % d == 0)
check("mollified coefficients: c(1)=1, c(n)=0 for 1<n<=X", cnum[1] == 1 and all(cnum[n] == 0 for n in range(2, Xm + 1)))


def MX(s):
    return sum(mu_tab[d] * chi4[d % 4] * mp.power(d, -s) for d in range(1, Xm + 1))


def lhs(s0):
    return sum(cnum[n] * chi4[n % 4] * mp.power(n, -s0) * mp.e ** (-mp.mpf(n) / Ym) for n in range(1, 2001))


def rhs_line(s0, c=-0.25):
    f = lambda v: mp.dirichlet(s0 + mp.mpc(c, v), chi4) * MX(s0 + mp.mpc(c, v)) * mp.gamma(mp.mpc(c, v)) * mp.power(Ym, mp.mpc(c, v))
    return mp.quad(f, [-30, -10, -3, 0, 3, 10, 30]) / (2 * mp.pi)


for s0, lab in [(rho1, "at the zero rho1 (residue L M = 0)"), (mp.mpc(0.8, 3.0), "at s0 = 0.8+3i (residue L M added)")]:
    L0 = lhs(s0)
    R0 = rhs_line(s0) + mp.dirichlet(s0, chi4) * MX(s0)
    check("detection identity sum c chi n^{-s} e^{-n/Y} = L M + line Re w = -1/4, " + lab, abs(L0 - R0) < 1e-8 * max(1, abs(L0)),
          "lhs %s, rhs %s" % (mp.nstr(L0, 10), mp.nstr(R0, 10)))
Y1m = Ym * math.log(Ym) ** 2
tail_num = sum(abs(cnum[n]) * n ** -0.5 * math.exp(-n / Ym) for n in range(int(Y1m) + 1, 2001))
check("numeric tail sum_{n>Y log^2 Y} |c(n)| n^{-1/2} e^{-n/Y} < 1e-3 (X=10, Y=30)", tail_num < 1e-3, "%.2e" % tail_num)

# partial summation inequality on random data
rng = random.Random(12345)
okp = True
worst = 0.0
for trial in range(300):
    U = rng.choice([20, 50, 100])
    sig = rng.uniform(0.5, 0.95)
    beta = rng.uniform(sig, 1.0)
    gam = rng.uniform(-50, 50)
    b = {n: complex(rng.gauss(0, 1), rng.gauss(0, 1)) for n in range(U + 1, 2 * U + 1)}
    piece = sum(b[n] * n ** complex(-beta, -gam) for n in b)
    S = 0j
    Svals = []
    for n in range(U + 1, 2 * U + 1):
        S += b[n] * n ** complex(-sig, -gam)
        Svals.append((n, S))
    S2U = Svals[-1][1]
    integ1 = 0.0
    integ2 = 0.0
    for idx, (n, Sn) in enumerate(Svals):    # S(u) = S_n on [n, n+1)
        hi = n + 1 if n < 2 * U else 2 * U
        w = math.log(hi / n) if hi > n else 0.0
        integ1 += abs(Sn) * w
        integ2 += abs(Sn) ** 2 * w
    bound1 = abs(S2U) + (beta - sig) * integ1
    bound2 = 2 * abs(S2U) ** 2 + math.log(2) / 2 * integ2
    supS = max(abs(Sn) for _, Sn in Svals)
    okp &= abs(piece) <= bound1 * (1 + 1e-12) and abs(piece) ** 2 <= bound2 * (1 + 1e-12) and abs(piece) ** 2 <= cps * supS ** 2 * (1 + 1e-12)
    worst = max(worst, abs(piece) / bound1)
check("partial summation: |sum b n^{-beta-i gamma}| <= |S(2U)| + (beta-sigma) int |S| du/u, and squared forms (300 random)", okp,
      "max ratio %.3f" % worst)

# Gallagher with multiplicity
okg = True
for trial in range(20):
    a = np.array([complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(40)])
    ln = np.log(np.arange(1, 41))
    centers = [rng.uniform(-8, 8) for _ in range(6)]
    pts = [c + rng.uniform(-0.15, 0.15) for c in centers for _ in range(rng.randint(1, 4))]
    tg = np.linspace(-10, 10, 40001)
    F = (a[None, :] * np.exp(-1j * np.outer(tg, ln))).sum(axis=1)
    dF = (a[None, :] * (-1j * ln)[None, :] * np.exp(-1j * np.outer(tg, ln))).sum(axis=1)
    G = np.abs(F) ** 2 + 2 * np.abs(F) * np.abs(dF)
    dt = tg[1] - tg[0]
    Zmax = max(sum(1 for p in pts if abs(p - t) <= 0.5) for t in np.linspace(-9, 9, 3601))
    lhs_g = sum(abs((a * np.exp(-1j * p * ln)).sum()) ** 2 for p in pts)
    rhs_g = Zmax * G[np.abs(tg) <= 9.0].sum() * dt
    okg &= lhs_g <= rhs_g
check("Gallagher with multiplicity: sum_rho |F(gamma)|^2 <= Z int (|F|^2 + 2|F||F'|) (20 random clustered sets)", okg)


# Lemma 5 on the coefficient vectors of the proof (modulus 2+i, chi_1 in the coefficients)
def gauss_factor(a, b):
    """prime-ideal factorisation of (a+bi) as list of (prime ideal key, exponent, norm)."""
    n = a * a + b * b
    fac = []
    m = n
    p = 2
    primes = []
    while p * p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            primes.append((p, e))
        p += 1
    if m > 1:
        primes.append((m, 1))
    for p, e in primes:
        if p == 2:
            fac.append((("2",), e, 2))
        elif p % 4 == 3:
            fac.append((("i", p), e // 2, p * p))
        else:
            # find pi = x + yi with x^2+y^2 = p, x > y > 0 ; count powers of pi and conj(pi) dividing a+bi
            x = next(xx for xx in range(1, math.isqrt(p) + 1) if math.isqrt(p - xx * xx) ** 2 == p - xx * xx)
            y = math.isqrt(p - x * x)
            for (u, v), lab in [((x, y), "+"), ((x, -y), "-")]:
                cnt = 0
                za, zb = a, b
                while True:
                    # (za+zb i)/(u+vi) = (za+zb i)(u-vi)/p
                    ra, rb = za * u + zb * v, zb * u - za * v
                    if ra % p == 0 and rb % p == 0:
                        za, zb = ra // p, rb // p
                        cnt += 1
                    else:
                        break
                if cnt:
                    fac.append((("s", p, lab), cnt, p))
    return fac


def mu_ideal(fac):
    if any(e > 1 for _, e, _ in fac):
        return 0
    return (-1) ** len(fac)


def c_moll(fac, X):
    """sum_{d | a, N d <= X} mu(d) over squarefree divisors."""
    ps = [nm for _, e, nm in fac]
    tot = 0
    for mask in range(1 << len(ps)):
        nd = 1
        sgn = 1
        for i, nm in enumerate(ps):
            if mask >> i & 1:
                nd *= nm
                sgn = -sgn
        if nd <= X:
            tot += sgn
    return tot


def lemma5_ratio(coef, ids, M, T, Nsup):
    th = np.array([math.atan2(b, a) for (a, b, n) in ids])
    lnN = np.log(np.array([n for (a, b, n) in ids], dtype=float))
    c = np.array(coef)
    dphi = 4 * (th[:, None] - th[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        DM = np.sin((M + 0.5) * dphi) / np.sin(dphi / 2)
    DM[np.abs(np.sin(dphi / 2)) < 1e-12] = 2 * M + 1
    dl = lnN[:, None] - lnN[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        KT = 2 * np.sin(T * dl) / dl
    KT[np.abs(dl) < 1e-14] = 2 * T
    lhs_v = np.real(np.conj(c) @ (DM * KT) @ c)
    cont = np.array([math.gcd(a, b) for (a, b, n) in ids])
    Nn = np.array([n for (a, b, n) in ids], dtype=float)
    D = cont * np.sqrt(Nsup / Nn)
    rhs_v = math.e ** 2 * math.sqrt(math.pi) * (math.sqrt(math.pi) * M + math.pi * Nsup / 2) * np.sum(np.abs(c) ** 2 * (T + math.sqrt(math.pi) * D))
    return lhs_v / rhs_v


Nsup = 400
idsL = ideal_gens(Nsup)
facs = {(a, b): gauss_factor(a, b) for (a, b, n) in idsL}
okf = all(np.prod([nm ** e for _, e, nm in facs[(a, b)]]) == n for (a, b, n) in idsL)
check("Gaussian factorisation consistent with norms (N a <= 400)", okf)
k_chi = 1
xi_chi = {(a, b): chi_mod_2pi(a, b, 1) * cmath.exp(1j * k_chi * math.atan2(b, a)) for (a, b, n) in idsL}
X5, Y5, sig5 = 10, 300.0, 0.75
worst = 0.0
for (U, u) in [(100, 180), (200, 400), (50, 100)]:
    sub = [(a, b, n) for (a, b, n) in idsL if U < n <= u]
    coefI = [c_moll(facs[(a, b)], X5) * xi_chi[(a, b)] * math.exp(-n / Y5) * n ** (-sig5) for (a, b, n) in sub]
    coefD = [cf * math.log(n) for cf, (a, b, n) in zip(coefI, sub)]
    for M in [1, 10, 100]:
        for T in [1.0, 5.0]:
            for cf in (coefI, coefD):
                if np.sum(np.abs(cf)) == 0:
                    continue
                worst = max(worst, lemma5_ratio(cf, sub, M, T, 2 * U))
check("Lemma 5 holds for the type-I vectors c(a) xi_chi(a) e^{-Na/Y} Na^{-sigma} on (U,u] and their log-derivatives", worst <= 1,
      "max LHS/bound %.4f" % worst)
worst = 0.0
for Np in [50, 100, 200, 400]:
    sub = [(a, b, n) for (a, b, n) in idsL if Np / 2 < n <= Np]
    coefM = [mu_ideal(facs[(a, b)]) * xi_chi[(a, b)] * n ** -0.5 for (a, b, n) in sub]
    for M in [1, 10, 100]:
        for T in [1.0, 5.0, 20.0]:
            worst = max(worst, lemma5_ratio(coefM, sub, M, T, Np))
check("Lemma 5 holds for the mollifier pieces mu(a) xi_chi(a) Na^{-1/2} on (N/2, N]", worst <= 1, "max LHS/bound %.4f" % worst)

# ---------------------------------------------------------------- E. ideal arithmetic
hdr("E. ideal arithmetic of Z[i]")
Nmax_e = 3000
ids_e = ideal_gens(Nmax_e)
r0 = [0] * (Nmax_e + 1)
for (a, b, n) in ids_e:
    r0[n] += 1
okr = all(r0[n] == sum(chi4[d % 4] for d in range(1, n + 1) if n % d == 0) for n in range(1, 1201))
check("#ideals of norm n = sum_{d|n} chi_-4(d) (n <= 1200)", okr)
dsum = [0] * (Nmax_e + 1)
for (a, b, n) in ids_e:
    fac = gauss_factor(a, b)
    dsum[n] += int(np.prod([e + 1 for _, e, _ in fac])) if fac else 1


def d4(n):
    # number of ordered factorisations into 4 factors
    tot = 1
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            tot *= math.comb(e + 3, 3)
        p += 1
    if m > 1:
        tot *= 4
    return tot


conv = [0] * (Nmax_e + 1)
for d in range(1, Nmax_e + 1):
    if r0[d]:
        for e in range(1, Nmax_e // d + 1):
            conv[d * e] += r0[d] * r0[e]
check("sum_{Na=n} d(a) = (r0*r0)(n) <= d_4(n) for n <= 3000", all(dsum[n] == conv[n] <= d4(n) for n in range(1, Nmax_e + 1)))
okd4 = True
acc = 0
for u in range(1, 100001):
    acc += d4(u)
    if u in (10, 100, 1000, 10000, 100000) or u % 9973 == 0:
        okd4 &= acc <= u * (1 + math.log(u)) ** 3
check("sum_{n<=u} d_4(n) <= u (1+log u)^3 (u <= 1e5, sampled)", okd4)
LamI = [0.0] * (Nmax_e + 1)
for (a, b, n) in ids_e:
    fac = gauss_factor(a, b)
    if len(fac) == 1:
        LamI[n] += math.log(fac[0][2])


def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0


check("sum_{Na=n} Lambda(a) = Lambda(n)(1 + chi_-4(n)) (n <= 3000): psi_K = psi + psi(chi_-4)",
      all(abs(LamI[n] - vonmangoldt(n) * (1 + chi4[n % 4])) < 1e-9 for n in range(2, Nmax_e + 1)))
lhs_c = sum(math.gcd(a, b) / n for (a, b, n) in ids_e)
rhs_c = (1 + math.log(math.sqrt(Nmax_e))) * sum(1 / n for (a, b, n) in ids_e)
check("sum_{Na<=X} cont(a)/Na <= (1+log sqrt X) sum 1/Na (X=3000)", lhs_c <= rhs_c, "%.2f <= %.2f" % (lhs_c, rhs_c))

# ---------------------------------------------------------------- F. sanity: primes in sectors
hdr("F. sanity only: Gaussian primes in classes and sectors, N z <= 1e7")
x = 10 ** 7
XM = int(1.04 * x)
sieve = np.ones(XM + 1, dtype=bool)
sieve[:2] = False
for p in range(2, math.isqrt(XM) + 1):
    if sieve[p]:
        sieve[p * p::p] = False
ppbase = np.zeros(XM + 1, dtype=np.int32)
pr = np.nonzero(sieve)[0]
ppbase[pr] = pr
for p in pr[pr <= math.isqrt(XM)]:
    q = int(p) * int(p)
    while q <= XM:
        ppbase[q] = p
        q *= int(p)
del sieve
R = math.isqrt(XM)
A_l, B_l, L_l = [], [], []
bb = np.arange(-R, R + 1, dtype=np.int64)
for a in range(-R, R + 1):
    n = a * a + bb * bb
    msk = n <= XM
    msk &= n > 0
    nn = n[msk]
    bsel = bb[msk]
    base = ppbase[nn]
    good = base > 0
    nn, bsel, base = nn[good], bsel[good], base[good]
    lam = np.zeros(len(nn))
    p1 = base % 4 == 1
    g = np.gcd(np.abs(bsel), abs(a))
    ok1 = p1 & (g % base != 0)
    lam[ok1] = np.log(base[ok1])
    lam[base == 2] = math.log(2)
    p3 = base % 4 == 3
    lam[p3] = 2 * np.log(base[p3])
    keep = lam > 0
    A_l.append(np.full(int(keep.sum()), a, dtype=np.int64))
    B_l.append(bsel[keep])
    L_l.append(lam[keep])
A_e = np.concatenate(A_l)
B_e = np.concatenate(B_l)
L_e = np.concatenate(L_l)
del ppbase, A_l, B_l, L_l
N_e = A_e * A_e + B_e * B_e
TH_e = np.arctan2(B_e, A_e)
print("elements with Lambda > 0 and N z <= 1.04e7: %d" % len(A_e))
inx = N_e <= x
psiK = L_e[inx].sum() / 4
check("prime ideal theorem sanity: psi_K(1e7)/1e7 within 5e-3 of 1", abs(psiK / x - 1) < 5e-3, "psi_K/x = %.6f" % (psiK / x))
r2pi = (A_e - 2 * B_e) % 5
r2mi = (A_e + 2 * B_e) % 5
anchors = [0.31, 0.77, 1.13, 2.40, -1.95, 3.00]
widths = [("x^-1/4", x ** -0.25), ("2x^-1/4", 2 * x ** -0.25), ("4x^-1/4", 4 * x ** -0.25), ("x^-0.27", x ** -0.27)]
zs = []
ratios_full = []
Lx = L_e[inx]
THx = TH_e[inx]
classes = [("2+i", r2pi[inx], None, 4)] + [("5", r2pi[inx], r2mi[inx], 16)]
for phi in anchors:
    for wl, om in widths:
        for orient in (+1, -1):
            lo = phi if orient > 0 else phi - om
            d = np.mod(THx - lo, 2 * np.pi)
            ins = (d > 0) & (d <= om)
            ratios_full.append(Lx[ins].sum() / ((2 / math.pi) * om * x))
            for name, c1, c2, Phi in classes:
                if c2 is None:
                    for r in range(1, 5):
                        sel = ins & (c1 == r)
                        S = Lx[sel].sum()
                        E = (2 / math.pi) * om * x / Phi
                        zs.append((S - E) / math.sqrt((Lx[sel] ** 2).sum()))
                else:
                    for r in range(1, 5):
                        for r2 in range(1, 5):
                            sel = ins & (c1 == r) & (c2 == r2)
                            S = Lx[sel].sum()
                            E = (2 / math.pi) * om * x / Phi
                            zs.append((S - E) / math.sqrt((Lx[sel] ** 2).sum()))
zs = np.array(zs)
rf = np.array(ratios_full)
print("full sectors (no class): Lambda-sum / ((2/pi) omega x): min %.4f max %.4f mean %.4f (%d sectors)" % (rf.min(), rf.max(), rf.mean(), len(rf)))
print("sector-class pairs: %d; z-scores mean %+.3f rms %.3f max|z| %.2f" % (len(zs), zs.mean(), math.sqrt((zs ** 2).mean()), np.abs(zs).max()))
check("main term sanity: full-sector ratios within 3%% of 1", np.all(np.abs(rf - 1) < 0.03))
check("main term sanity: class z-scores |z| < 4.5 (960 pairs)", np.abs(zs).max() < 4.5)

# sizes of Psi(xi) with the log-Gaussian weight (sanity only)
H = math.log(x) ** 2
lu = np.log(N_e / x) * H
V = np.where(lu < -8, 1.0, np.where(lu > 8, 0.0, 0.5))
mid = np.abs(lu) <= 8
V[mid] = np.array([0.5 * math.erfc(t / math.sqrt(2)) for t in lu[mid]])
u = np.exp(1j * TH_e)
w = L_e * V / 4
res = {}
for (mod, j) in [(1, 0), (5, 0), (5, 1), (5, 2), (5, 3)]:
    if mod == 1:
        base = w.astype(complex)
    else:
        chv = np.zeros(len(A_e), dtype=complex)
        for r in range(1, 5):
            chv[r2pi == r] = 1j ** (j * LOG5[r])
        base = w * chv
    worst = 0.0
    zero_ok = True
    zk = base * u ** (-200)
    step = u
    for k in range(-200, 201):
        val = zk.sum()
        if (k - j) % 4 == 0:
            if not (mod in (1, 5) and j == 0 and k == 0):
                worst = max(worst, abs(val))
        else:
            zero_ok &= abs(val) < 1e-6 * x
        zk = zk * step
    res[(mod, j)] = worst
    check("unit sum: element sum vanishes for k != k_chi mod 4 (modulus %s, chi_%d)" % ("1" if mod == 1 else "2+i", j), zero_ok)
    print("   modulus %s chi_%d: max over admissible 0<|k|<=200 (and k=0 if chi nontrivial) of |Psi| = %.1f = %.2e x = %.2f sqrt(x)"
          % ("1" if mod == 1 else "2+i", j, worst, worst / x, worst / math.sqrt(x)))
check("sanity: max |Psi(xi)| over the nontrivial family |k|<=200 is <= 1e-2 x", max(res.values()) <= 1e-2 * x)
del A_e, B_e, L_e, N_e, TH_e, V, u, w

# ---------------------------------------------------------------- G. illustration
hdr("G. illustration (heuristic; placeholder constants c = delta' = 0.1, B = 30): VK vs classical")
print("   log10 of the zero-sum bound divided by x: x^{-delta' eta} (log x)^B")
def g_vk(lx):
    return (-0.1 * 0.1 * lx ** (1 / 3) * math.log(lx) ** (-1 / 3) + 30 * math.log(lx)) / math.log(10)


def g_cl(lx):
    return (-0.1 * 0.1 + 30 * math.log(lx)) / math.log(10)


for lx in [1e3, 1e9, 1e15, 1e18, 1e21, 1e24]:
    print("   log x = %.0e: VK %+.3e   classical %+.3e" % (lx, g_vk(lx), g_cl(lx)))
check("illustration: eta = c(log x)^{-2/3}(loglog x)^{-1/3}: bound/x -> 0; eta = c/log x: bound/x -> infinity",
      g_vk(1e24) < g_vk(1e21) < 0 and g_cl(1e24) > g_cl(1e21) > 0)

# ---------------------------------------------------------------- H. repaired inequalities
hdr("H. inequalities repaired after the referee report")
T_H = time.time()


# H1. consequence (a) of Lemma PF --------------------------------------------------------------------------------
# The proof uses only that every window {|gamma - T| <= 1} holds at most the stated number of zeros, so the actual
# window counts of a synthetic ordinate set may stand in for ZC(T). Both signs eps = +1, -1 (reflect the ordinates).
def ZCf(Nf, k, T):
    return 1.25 * math.log(Nf) + 2.5 * math.log(abs(k) + abs(T) + 2) + 3.35


def wcount(gs, T):
    return int(np.count_nonzero(np.abs(gs - T) <= 1.0))


def pf_case(gs, n, T=None, rng_=None):
    """choose T in [n, n+1] at distance >= 1/(2 Z_n + 2) from every ordinate (Z_n = window count at n + 1/2);
    return (T found?, S = sum_{|gamma-T|<1} 1/|T-gamma|, new bound, old bound, some near zero outside [n-1/2, n+3/2]?)"""
    Zn = wcount(gs, n + 0.5)
    r = 1.0 / (2 * Zn + 2)
    if T is None:
        grid = np.linspace(n, n + 1, 2001)
        dmin = np.min(np.abs(grid[:, None] - gs[None, :]), axis=1) if len(gs) else np.full(len(grid), 9.0)
        good = np.nonzero(dmin >= r)[0]
        if len(good) == 0:
            return False, 0.0, 0.0, 0.0, False
        T = float(grid[good[rng_.randrange(len(good))]])
    near = gs[np.abs(gs - T) < 1]
    S = float(np.sum(1.0 / np.abs(T - near)))
    new = Zn * (2 * Zn + 2) + 2 * (wcount(gs, n - 0.5) + wcount(gs, n + 1.5))
    old = Zn * (2 * Zn + 2)
    outside = bool(np.any((near < n - 0.5) | (near > n + 1.5)))
    return True, S, new, old, outside


rngH = random.Random(20261126)
ok_found = ok_new = True
n_cases = n_outside = n_oldfail = 0
worst_ratio = 0.0
for trial in range(300):
    cents = [rngH.uniform(-14, 14) for _ in range(rngH.randint(3, 30))]
    gs0 = np.array([c + rngH.uniform(-0.3, 0.3) for c in cents for _ in range(rngH.randint(1, 3))])
    for eps_sign in (1, -1):
        gs = eps_sign * gs0
        for n in range(1, 12):
            found, S, new, old, outside = pf_case(gs, n, rng_=rngH)
            ok_found &= found
            if not found:
                continue
            n_cases += 1
            ok_new &= S <= new * (1 + 1e-12)
            if new > 0:
                worst_ratio = max(worst_ratio, S / new)
            n_outside += outside
            n_oldfail += S > old
check("Lemma PF (a): T_n in [n, n+1] at distance >= 1/(2Z_n+2) from all ordinates exists (both signs, %d cases)" % n_cases, ok_found)
check("Lemma PF (a): sum_{|gamma-T_n|<1} 1/|T_n-gamma| <= Z_n(2Z_n+2) + 2(ZC(n-1/2) + ZC(n+3/2)) (window counts as ZC)",
      ok_new, "max ratio %.3f; cases with a zero at |gamma-T_n|<1 outside [n-1/2, n+3/2]: %d of %d; old bound exceeded: %d"
      % (worst_ratio, n_outside, n_cases, n_oldfail))
# negative control: one zero at n - 0.6 (outside [n-1/2, n+3/2]), T_n = n: the old bound Z_n(2Z_n+2) = 0 fails
_, S_nc, new_nc, old_nc, out_nc = pf_case(np.array([3 - 0.6]), 3, T=3.0)
check("Lemma PF (a), negative control {n - 0.6}, T_n = n: old bound fails, new bound holds", S_nc > old_nc and S_nc <= new_nc and out_nc,
      "S = %.4f, old %.1f, new %.1f" % (S_nc, old_nc, new_nc))
# second line of (a): with ZC(.) <= ZC(n+3/2) <= 4.6 l, the first line is <= 60 l^2 + C_PF l, l = log(Nf(|k|+n+4))
ok_l = True
worst_l = 0.0
for Nf in [1, 2, 5, 13, 100, 10 ** 4, 10 ** 8]:
    for k in [0, 1, 3, 10, 100, 10 ** 4, 10 ** 6]:
        for n in [1, 2, 3, 5, 10, 100, 10 ** 4, 10 ** 6]:
            ell = math.log(Nf * (abs(k) + n + 4))
            Zn = ZCf(Nf, k, n + 0.5)
            first = Zn * (2 * Zn + 2) + 2 * (ZCf(Nf, k, n - 0.5) + ZCf(Nf, k, n + 1.5))
            ok_l &= ell > 1.6 and ZCf(Nf, k, n + 1.5) <= 4.6 * ell and first <= 60 * ell ** 2
            worst_l = max(worst_l, first / ell ** 2)
check("Lemma PF (a), second line: ZC(n+3/2) <= 4.6 l and Z_n(2Z_n+2) + 2(ZC(n-1/2)+ZC(n+3/2)) <= 60 l^2 (grid Nf, k, n)", ok_l,
      "max first-line/l^2 = %.2f" % worst_l)

# H2. (T1) and (T2) of section 4.3 on a complete range -----------------------------------------------------------
NT = 10 ** 6


def conv(f, g):
    """Dirichlet convolution (f*g)(n), n <= NT."""
    out = np.zeros(NT + 1, dtype=np.int64)
    for m in range(1, NT + 1):
        fm = int(f[m])
        if fm:
            out[m::m] += fm * g[1:NT // m + 1]
    return out


one = np.ones(NT + 1, dtype=np.int64)
one[0] = 0
dj = {1: one}
for j in range(2, 5):
    dj[j] = conv(dj[j - 1], one)
uu = np.arange(1, NT + 1, dtype=np.float64)
lg = np.log(uu)
harm = np.cumsum(1.0 / uu)
check("(T1) sum_{m<=u} 1/m <= 1 + log u for every integer u <= 1e6", bool(np.all(harm <= 1 + lg + 1e-12)),
      "max excess %.3e" % float(np.max(harm - 1 - lg)))
ok_t1 = True
info = []
for j in range(1, 5):
    cs = np.cumsum(dj[j][1:]).astype(np.float64)
    rhs = uu * (1 + lg) ** (j - 1)
    ok_t1 &= bool(np.all(cs <= rhs * (1 + 1e-12)))
    info.append("j=%d max ratio %.3f" % (j, float(np.max(cs / rhs))))
check("(T1) induction: sum_{n<=u} d_j(n) <= u(1+log u)^(j-1), j = 1..4, every integer u <= 1e6 (so every real u in [1, 1e6+1))",
      ok_t1, "; ".join(info))
ok_d4 = all(int(dj[4][n]) == d4(n) for n in range(1, 3001)) and all(int(dj[4][n]) == d4(n) for n in rngH.sample(range(3001, NT + 1), 300))
check("d_4 = 1*1*1*1 by convolution agrees with the factorisation formula prod C(e+3,3) (n <= 3000 and 300 random n <= 1e6)", ok_d4)
dd = conv(dj[2], dj[2])
check("d * d = d_4 for n <= 1e6", bool(np.array_equal(dd, dj[4])))
del dd
# r0 by counting lattice points: ideals of norm n <-> a + bi with a >= 1, b >= 0
Rr = math.isqrt(NT)
norms = []
for a in range(1, Rr + 1):
    bmax = math.isqrt(NT - a * a)
    norms.append(a * a + np.arange(0, bmax + 1, dtype=np.int64) ** 2)
r0b = np.bincount(np.concatenate(norms), minlength=NT + 1).astype(np.int64)
del norms
chiv = np.zeros(NT + 1, dtype=np.int64)
chiv[1::4] = 1
chiv[3::4] = -1
check("r0 = 1 * chi_-4 (lattice count) for n <= 1e6", bool(np.array_equal(r0b, conv(chiv, one))))
check("(T2) r0(n) <= d(n) for every n <= 1e6", bool(np.all(r0b[1:] <= dj[2][1:])))
ok_pp = True
for p in range(2, 1001):
    if all(p % q for q in range(2, math.isqrt(p) + 1)):
        e, pe = 1, p
        while pe <= NT:
            want = 1 if p == 2 else (e + 1 if p % 4 == 1 else (1 if e % 2 == 0 else 0))
            ok_pp &= int(r0b[pe]) == want
            e += 1
            pe *= p
ok_mult = all(int(r0b[m * n]) == int(r0b[m]) * int(r0b[n]) for m, n in
              ((rngH.randint(1, 1000), rngH.randint(1, 1000)) for _ in range(3000)) if math.gcd(m, n) == 1)
check("(T2) r0(2^e) = 1, r0(p^e) = e+1 (p = 1 mod 4), 1 or 0 (p = 3 mod 4, e even/odd), p <= 1000; r0 multiplicative (random pairs)",
      ok_pp and ok_mult)
rr = conv(r0b, r0b)
check("(T2) (r0*r0)(n) <= d_4(n) for every n <= 1e6", bool(np.all(rr[1:] <= dj[4][1:])))
csr = np.cumsum(rr[1:]).astype(np.float64)
check("(T2) sum_{Na<=u} d(a) = sum_{n<=u} (r0*r0)(n) <= u(1+log u)^3 for every integer u <= 1e6", bool(np.all(csr <= uu * (1 + lg) ** 3)),
      "max ratio %.3f" % float(np.max(csr / (uu * (1 + lg) ** 3))))
# sum_{Na=n} d(a) by factorising every ideal of norm <= 22000 (independent of the convolution identity)
NMAXT = 22000
dsum_big = np.zeros(NMAXT + 1, dtype=np.int64)
for (a, b, n) in ideal_gens(NMAXT):
    fac = gauss_factor(a, b)
    dsum_big[n] += int(np.prod([e + 1 for _, e, _ in fac])) if fac else 1
check("(T2) sum_{Na=n} d(a) (ideal factorisation) = (r0*r0)(n) for n <= 22000", bool(np.array_equal(dsum_big[1:], rr[1:NMAXT + 1])))
del dj, rr, csr, r0b, chiv, one, uu, lg, harm


# H3. (T3): the dyadic tail bound for every Y >= 100 ----------------------------------------------------------------
def log_tj(j, L):
    return (1 + j / 2) * math.log(2) + np.log(L) + L / 2 + 3 * np.log(1 + (j + 1) * math.log(2) + L + 2 * np.log(L)) - 2.0 ** j * L ** 2


def dlog_tj(j, L):
    return 0.5 + 1 / L + 3 * (1 + 2 / L) / (1 + (j + 1) * math.log(2) + L + 2 * np.log(L)) - 2.0 ** (j + 1) * L


Lg = np.linspace(4.0, 200.0, 200001)
ok_mono = all(bool(np.all(dlog_tj(j, Lg) <= 0.5 + 0.25 + 0.9 - 2.0 ** (j + 1) * Lg)) and bool(np.all(dlog_tj(j, Lg) < 0))
              and bool(np.all(np.diff(log_tj(j, Lg)) < 0)) for j in range(0, 25))
check("(T3) t_j(L) decreasing on L >= 4 (derivative formula <= 1.65 - 2^(j+1) L < 0, and on a grid of 2e5 points), j = 0..24", ok_mono)
Lh = Lg[Lg >= math.log(100)]
rat = max(float(np.max(log_tj(j + 1, Lh) - log_tj(j, Lh))) for j in range(0, 24))
rb = math.sqrt(2) * (1 + math.log(2) / (1 + math.log(100))) ** 3 * math.exp(-math.log(100) ** 2)
check("(T3) t_(j+1)/t_j <= sqrt2 (1 + log2/(1+L))^3 e^(-L^2) <= 1.3e-9 for L >= log 100", math.exp(rat) <= rb * (1 + 1e-9) and rb <= 1.3e-9,
      "max ratio %.3e, bound %.3e" % (math.exp(rat), rb))
t0_100 = math.exp(float(log_tj(0, np.array([math.log(100)]))[0]))
L100 = math.log(100)
check("(T3) t_0(log 100) = 2*10*L(1 + log2 + L + 2 log L)^3 e^(-L^2) = 4.642e-5 (base 9.3527, L^2 = 21.2076)",
      abs(t0_100 - 4.642e-5) < 5e-9 and abs(1 + math.log(2) + L100 + 2 * math.log(L100) - 9.3527) < 1e-4 and abs(L100 ** 2 - 21.2076) < 1e-4,
      "t_0 = %.6e" % t0_100)


def dyadic_bound(Y):
    L = math.log(Y)
    return sum(math.exp(float(log_tj(j, np.array([L]))[0])) for j in range(0, 40))


Ygrid = np.concatenate([np.linspace(100, 1000, 901), np.logspace(3, 300, 600)])
bvals = np.array([dyadic_bound(Y) for Y in Ygrid])
check("(T3) sum_j t_j(log Y) <= 4.7e-5 for Y in [100, 1e300] (1501 points) and non-increasing in Y",
      bool(np.all(bvals <= 4.7e-5)) and bool(np.all(np.diff(bvals) <= 0)), "max %.4e at Y = %.0f" % (bvals.max(), Ygrid[bvals.argmax()]))
check("(T3) the closed form equals the dyadic sum tail_bound of section B (Y = 100, 1e3, 1e6)",
      all(abs(dyadic_bound(Y) - tail_bound(Y)) <= 1e-12 * tail_bound(Y) + 1e-300 for Y in [100.0, 1e3, 1e6]))
ok_true = True
info = []
for Y in [100.0, 150.0, 200.0, 300.0]:
    Y1 = Y * math.log(Y) ** 2
    nn = np.arange(int(Y1) + 1, NMAXT + 1)
    direct = float(np.sum(dsum_big[nn] * nn ** -0.5 * np.exp(-nn / Y)))
    rem = sum(2 * (2 ** j * NMAXT) ** 0.5 * (1 + math.log(2 ** (j + 1) * NMAXT)) ** 3 * math.exp(-(2 ** j) * NMAXT / Y) for j in range(0, 60))
    ok_true &= direct + rem <= dyadic_bound(Y) <= 1e-4
    info.append("Y=%g: tail %.2e (+rem %.1e) <= %.2e" % (Y, direct, rem, dyadic_bound(Y)))
check("(T3) true tail sum_{Na>Y1} d(a) Na^(-1/2) e^(-Na/Y) (ideals to 22000, remainder by (T2)) <= dyadic bound <= 1e-4",
      ok_true, "; ".join(info))
del dsum_big

# H4. the bound (Gamma) of Lemma DI, through the shift; (C3)(i), (iii) ----------------------------------------------
mp.mp.dps = 20
vg = np.concatenate([np.linspace(0, 0.95, 20), np.linspace(1, 10, 91), np.linspace(11, 100, 90), [150.0, 300.0, 500.0]])
avals = [round(-0.5 + 0.05 * i, 10) for i in range(51)]
wG = wHi = wLo = wSh = 0.0
for a in avals:
    for v0 in vg:
        for v in (float(v0), -float(v0)):
            if a == 0 and v == 0:
                continue          # the pole w = 0
            w = mp.mpc(a, v)
            g = abs(mp.gamma(w))
            av = abs(v)
            dec = mp.e ** (-mp.pi * av / 2)
            if a >= 0.5:
                wHi = max(wHi, float(g / (mp.sqrt(2 * mp.pi) * mp.e ** (1 / (12 * mp.mpf(a))) * abs(w) ** (mp.mpf(a) - 0.5) * dec)))
            if av >= 1:
                wG = max(wG, float(g / (10 * mp.mpf(av) ** 1.5 * dec)))
                if a < 0.5:
                    wLo = max(wLo, float(g / (mp.mpf(5.34) * dec)))
            if a < 0.5:
                wSh = max(wSh, float(abs(mp.gamma(w) - mp.gamma(w + 1) / w) / g))
check("(Gamma) |Gamma(a+iv)| <= sqrt(2pi) e^(1/(12a)) |w|^(a-1/2) e^(-pi|v|/2), 1/2 <= a <= 2, all v on the grid", wHi <= 1,
      "max ratio %.4f" % wHi)
check("(Gamma) shifted case: |Gamma(a+iv)| <= 5.34 e^(-pi|v|/2) for -1/2 <= a < 1/2, |v| >= 1 (via Gamma(w) = Gamma(w+1)/w)", wLo <= 1,
      "max ratio %.4f" % wLo)
check("(Gamma) |Gamma(a+iv)| <= 10 |v|^(3/2) e^(-pi|v|/2) for -1/2 <= a <= 2, 1 <= |v| <= 500 (51 x 368 grid)", wG <= 1,
      "max ratio %.4f" % wG)
check("(C3)(i) shift Gamma(w) = Gamma(w+1)/w on -1/2 <= Re w < 1/2 (relative error)", wSh < 1e-15, "max %.1e" % wSh)
shape = float(abs(mp.gamma(mp.mpc(2, 1e4))) * mp.e ** (mp.pi * 1e4 / 2) / mp.mpf(1e4) ** 1.5)
nc_g = float(abs(mp.gamma(mp.mpc(2, 100))) * mp.e ** (mp.pi * 50))
check("(Gamma) shape: at Re w = 2, |Gamma| e^(pi|v|/2)/|v|^(3/2) -> sqrt(2 pi); negative control: 10 e^(-pi|v|/2) without |v|^(3/2) fails at 2+100i",
      abs(shape - math.sqrt(2 * math.pi)) < 1e-3 and nc_g > 10, "ratio at v=1e4: %.5f; |Gamma(2+100i)| e^(50 pi) = %.1f" % (shape, nc_g))
ok_cm = True
for y in [0.3, 1.7]:
    for c in [0.5, 2.0]:
        val = mp.quad(lambda v: mp.gamma(mp.mpc(c, v)) * mp.power(y, -mp.mpc(c, v)), [-80, -20, -5, 0, 5, 20, 80]) / (2 * mp.pi)
        ok_cm &= abs(val - mp.e ** (-y)) < 1e-12
check("(C3)(iii) Cahen-Mellin e^(-y) = (1/2 pi i) int_(c) Gamma(w) y^(-w) dw (y = 0.3, 1.7; c = 1/2, 2)", ok_cm)
check("(C3)(i) Gamma(1) = 1 and residue of Gamma at 0 equals 1", mp.gamma(1) == 1 and abs(mp.mpf("1e-12") * mp.gamma(mp.mpf("1e-12")) - 1) < 1e-11)

# H5. blocks in 1 - sigma (section 5): 2^(A(j+1)/log x) and the block sum -------------------------------------------
ok_bl = True
for kap in [0.01, 0.1, 1 / 6, 103 / 512, 0.2499]:
    A = 2 + 8 * kap
    for lx in [4.0, 5.0, 10.0, 50.0, 100.0, 1e3, 1e4, 1e6]:
        jmax = math.floor(lx / A)
        ok_bl &= 2 ** (A * (jmax + 1) / lx) <= 4
        for lam in [0.05, 0.9 / A, 0.999 / A]:
            if lam * lx > 700:
                continue
            M = math.ceil(math.exp(lam * lx))
            dp = 1 - lam * A
            for j in range(0, jmax + 1, max(1, jmax // 40)):
                lhs_b = A * (j + 1) / lx * math.log(M)
                rhs_b = A * (j + 1) / lx * math.log(2) + (1 - dp) * (j + 1)
                ok_bl &= lhs_b <= rhs_b + 1e-9
            ok_bl &= math.log(M) <= lx and lx ** 3 + 4 <= 2 * lx ** 3
check("blocks: M^(A(j+1)/log x) <= 2^(A(j+1)/log x) e^((1-delta')(j+1)), and 2^(A(j+1)/log x) <= 4 for j <= log x/A (grid)", ok_bl)
ok_gs = True
for dp in [1e-3, 0.01, 0.1, 0.5, 1.0]:
    for etaL in np.linspace(0, 300, 601):
        j0 = math.floor(etaL)
        ok_gs &= math.exp(-dp * j0) / (1 - math.exp(-dp)) <= (2 * math.e / dp) * math.exp(-dp * etaL) * (1 + 1e-12)
check("blocks: sum_{j >= floor(eta L)} e^(-delta' j) <= (2e/delta') e^(-delta' eta L) (grid delta', eta L)", ok_gs)
print("   section H: %.1f s" % (time.time() - T_H))

print("\n%d checks, %d failed  (%.1f s)" % (NCHK, NFAIL, time.time() - T_START))
