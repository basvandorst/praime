"""Rigorous ball-arithmetic tools for ONE element f = sum_j b_j g_j of the family (b_j exact rationals).
This is certs/certify.py rewritten as a class (same mathematics, same bounds, same constants), so that it can be
applied to f = gt - eta h and to h itself.  See certify.py / README.md (section 4) for the derivations:
  T_N(u) = 4 d(N) (pi N)^{-1/2} X^{1/2} [P0(X) K0(2X) + P1(X) K1(2X)],  X = pi N e^u,  f = sum_N T_N,
  |T_N^{(r)}| <= h(N) = 4.4 Pabs_r(1) X^M e^{-2X}, ratios h(N+1)/h(N) <= e^{1-2pi} < 0.0051 for N >= M;
  K0, K1 from kbessel.K01 (series / asymptotic with DLMF 10.40.10 remainder).
Methods: point_vals(u, rs), crude_sup(r, a, b), lower_bound_interval(a, b), positive_on(intervals),
tail_positive(U) (asymptotic argument on [U, oo)), cost_ball() (2 sum Lambda(n) n^{-1/2} f(log n), balls n <= 80
plus explicit tails), value0()."""
import time
from fractions import Fraction
import common  # noqa: F401
from flint import arb, arb_poly, fmpq, ctx
from core import basis_polys, divisor_count, von_mangoldt
from kbessel import K01

_BP = {}


def fq(x):
    return arb(fmpq(x.numerator, x.denominator))


def strip(c):
    c = list(c)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def a_coef(nu, k):
    num = Fraction(1)
    for i in range(1, k + 1):
        num *= (4 * nu * nu - (2 * i - 1) ** 2)
    f = 1
    for i in range(1, k + 1):
        f *= i
    return num / (f * 8 ** k)


class Cert:
    def __init__(self, b, PT=16):
        self.b = [Fraction(x) for x in b]
        d = len(b) - 1
        self.d, self.PT = d, PT
        RMAX = PT + 1
        key = (d, RMAX)
        if key not in _BP:
            _BP[key] = basis_polys(d, RMAX)
        BP = _BP[key]
        n = d + 1
        P = []
        for r in range(RMAX + 1):
            L0 = max(len(BP[j][r][0]) for j in range(n))
            L1 = max(len(BP[j][r][1]) for j in range(n))
            P0 = [sum((self.b[j] * BP[j][r][0][k] for j in range(n) if k < len(BP[j][r][0])), Fraction(0)) for k in range(L0)]
            P1 = [sum((self.b[j] * BP[j][r][1][k] for j in range(n) if k < len(BP[j][r][1])), Fraction(0)) for k in range(L1)]
            P.append((strip(P0), strip(P1)))
        self.P = P
        self.AP = [(arb_poly([fq(x) for x in P0]), arb_poly([fq(x) for x in P1])) for (P0, P1) in P]
        self.APabs = [(arb_poly([fq(abs(x)) for x in P0]), arb_poly([fq(abs(x)) for x in P1])) for (P0, P1) in P]
        self.Mdeg = [max(len(P0), len(P1)) - 1 for (P0, P1) in P]
        self.Pabs1 = [fq(sum(abs(x) for x in P0) + sum(abs(x) for x in P1)) for (P0, P1) in P]
        self.PI = arb.pi()
        self.RATIO = arb("0.0051")
        assert (self.RATIO - (1 - 2 * self.PI).exp()).lower() > 0
        self.FACT = [arb(1)]
        for k in range(1, RMAX + 2):
            self.FACT.append(self.FACT[-1] * k)

    # ---- identical to certify.py ----
    def tail_bound(self, r, a, bnd, N0):
        M = self.Mdeg[r]
        ea = fq(a).exp()
        eb = fq(bnd).exp()
        N1 = max(N0, M) + 1
        tot = arb(0)
        for N in range(N0 + 1, N1 + 1):
            h = arb("4.4") * self.Pabs1[r] * (self.PI * N * eb) ** M * (-2 * self.PI * N * ea).exp()
            tot += h / (1 - self.RATIO) if N == N1 else h
        return tot.upper()

    def point_vals(self, u, rs):
        PI = self.PI
        eu = fq(u).exp()
        out = {r: arb(0) for r in rs}
        N = 1
        while True:
            X = PI * N * eu
            k0, k1 = K01(2 * X)
            pref = 4 * divisor_count(N) / (PI * N).sqrt() * X.sqrt()
            for r in rs:
                out[r] += pref * (self.AP[r][0](X) * k0 + self.AP[r][1](X) * k1)
            if N == 1:
                T1 = abs(pref * (self.APabs[rs[0]][0](X) * k0 + self.APabs[rs[0]][1](X) * k1))
            tb = max(self.tail_bound(r, u, u, N) for r in rs)
            if tb < (T1 * arb("1e-70")).lower():
                break
            N += 1
        for r in rs:
            out[r] += arb(0, self.tail_bound(r, u, u, N))
        return out

    def crude_sup(self, r, a, bnd):
        PI = self.PI
        ea = fq(a).exp()
        eb = fq(bnd).exp()
        tot = arb(0)
        N = 1
        while True:
            Xa = PI * N * ea
            Xb = PI * N * eb
            k0, k1 = K01(2 * Xa)
            term = 4 * divisor_count(N) / (PI * N).sqrt() * Xb.sqrt() * (self.APabs[r][0](Xb) * k0 + self.APabs[r][1](Xb) * k1)
            tot += term
            if N == 1:
                T1 = term
            tb = self.tail_bound(r, a, bnd, N)
            if tb < (T1 * arb("1e-70")).lower():
                break
            N += 1
        return (tot + arb(self.tail_bound(r, a, bnd, N))).upper()

    def lower_bound_interval(self, a, bnd):
        PT = self.PT
        c = (a + bnd) / 2
        rad = (bnd - a) / 2
        ra = fq(rad)
        v = self.point_vals(c, list(range(PT + 1)))
        g0, g1, g2 = v[0], v[1], v[2]
        lo_q = min((g0 + g1 * ra + g2 * ra * ra / 2).lower(), (g0 - g1 * ra + g2 * ra * ra / 2).lower())
        if g2.lower() > 0:
            x = -g1 / g2
            if not (abs(x).lower() > ra.upper()):
                lo_q = min(lo_q, (g0 - g1 * g1 / (2 * g2)).lower())
        elif not (g2.upper() < 0):
            lo_q = min(lo_q, (g0 - abs(g1) * ra - abs(g2) * ra * ra / 2).lower())
        hi = arb(0)
        for k in range(3, PT + 1):
            hi += abs(v[k]) * ra ** k / self.FACT[k]
        rem = arb(self.crude_sup(PT + 1, a, bnd)) * ra ** (PT + 1) / self.FACT[PT + 1]
        lo = arb(lo_q) - hi - rem
        return lo, g0

    def positive_on(self, intervals, min_width=Fraction(1, 2 ** 70), log=print, t0=None):
        """prove f > 0 on every closed interval of the list (Fractions); bisection like certify.py.
        Returns dict(intervals=..., depth_max=..., worst=(rel, a, b))."""
        t0 = time.time() if t0 is None else t0
        stack = list(reversed(intervals))
        nint = 0
        worst = None
        depth_max = 0
        w0 = {}
        while stack:
            a, bnd = stack.pop()
            lo, gc = self.lower_bound_interval(a, bnd)
            nint += 1
            if lo.lower() > 0:
                rel = float((lo / gc).lower())
                if worst is None or rel < worst[0]:
                    worst = (rel, float(a), float(bnd))
                depth_max = max(depth_max, int(Fraction(1, 64) / (bnd - a)).bit_length() - 1)   # floor(log2((1/64)/width))
                if nint % 500 == 0:
                    log("    ... %d intervals, at u=%.6f, %.0fs" % (nint, float(bnd), time.time() - t0))
                continue
            if bnd - a < min_width:
                raise RuntimeError("FAILED to prove f > 0 near u = %s (f(c) = %s, lower bound %s)" % (float(a), gc, lo))
            m = (a + bnd) / 2
            stack.append((m, bnd))
            stack.append((a, m))
        return dict(intervals=nint, depth_max=depth_max, worst=worst)

    def tail_positive(self, U, ELL=6):
        """certify.py part (3): f > 0 on [U, oo) via the asymptotic expansion of K0, K1 (DLMF 10.40.10)."""
        PI = self.PI
        d = self.d
        P0, P1 = self.P[0]
        M = len(P0) - 1
        assert M == 4 + 2 * d and len(P1) - 1 == 3 + 2 * d, "top coefficient vanishes"
        XU = PI * fq(U).exp()
        y = 1 / XU
        theta = arb("0.5 +/- 0.5")
        coef = {}
        for (Pc, nu) in ((P0, 0), (P1, 1)):
            for k, pk in enumerate(Pc):
                if pk == 0:
                    continue
                for i in range(ELL + 1):
                    c = fq(pk) * fq(a_coef(nu, i)) / arb(2) ** i
                    if i == ELL:
                        c = c * theta
                    e = (M - k) + i
                    coef[e] = coef.get(e, arb(0)) + c
        c0 = coef.pop(0)
        rest = arb(0)
        for e, c in coef.items():
            rest += abs(c) * y ** e
        Flow = c0 - rest
        ratio = arb("2.2") * self.Pabs1[0] * arb(2) ** M * (-2 * XU).exp() / (Flow * (1 - arb(2) ** M * (-2 * XU).exp()))
        ok = c0.lower() > 0 and Flow.lower() > 0 and ratio.upper() < 1
        return ok, dict(tail_leading=c0.str(10), tail_Flow=Flow.str(10), tail_ratio=ratio.str(5))

    def value0(self):
        return self.point_vals(Fraction(0), [0])[0]

    def cost_ball(self):
        PI = self.PI
        AP, APabs, Mdeg, Pabs1, RATIO = self.AP, self.APabs, self.Mdeg, self.Pabs1, self.RATIO
        cost = arb(0)
        NMAXC = max(80, Mdeg[0] + 2)
        for m in range(2, NMAXC + 1):
            lam = von_mangoldt(m)
            if lam == 0:
                continue
            val = arb(0)
            N = 1
            while True:
                X = PI * N * m
                k0, k1 = K01(2 * X)
                pref = 4 * divisor_count(N) / (PI * N).sqrt() * X.sqrt()
                val += pref * (AP[0][0](X) * k0 + AP[0][1](X) * k1)
                if N == 1:
                    T1 = abs(pref * (APabs[0][0](X) * k0 + APabs[0][1](X) * k1))
                tb = arb(0)
                N1 = max(N, Mdeg[0]) + 1
                for NN in range(N + 1, N1 + 1):
                    h = arb("4.4") * Pabs1[0] * (PI * NN * m) ** Mdeg[0] * (-2 * PI * NN * m).exp()
                    tb += h / (1 - RATIO) if NN == N1 else h
                tb = tb.upper()
                if tb < (T1 * arb("1e-70")).lower():
                    break
                N += 1
            val += arb(0, tb)
            p_ = next(p for p in range(2, m + 1) if m % p == 0)
            cost += 2 * arb(p_).log() / arb(m).sqrt() * val
        Xn = PI * (NMAXC + 1)
        tailc = 2 * arb(NMAXC + 1).log() * arb("4.4") * Pabs1[0] * Xn ** Mdeg[0] * (-2 * Xn).exp() / (1 - RATIO) ** 3
        cost += arb(0, tailc.upper())
        return cost
