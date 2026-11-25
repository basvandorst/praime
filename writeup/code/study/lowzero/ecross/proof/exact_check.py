"""Exact identities (no stationary phase), Gaussian-edged box H = 1_[T1,T2] * g_Delta:
 (A) Prop 1:  sum_gamma H(gamma) e^{2 i theta(gamma)} = -(1/2pi) sum_n Lambda_E(n) n^{-1/2} I_n  (+ E_0, negative zeros: < 1e-6)
 (B) Voronoi: sum_n lambda(n) n^{-1/2} I_n = w sum_m lambda(m) m^{-1/2} Hhat(-ln m),  Hhat(-ln m) = int H(t) m^{it} dt
with I_n = int H(t) e^{i(2 theta(t) - t ln n)} dt by the trapezoid rule (Horner in z = e^{-i h ln n})."""
import os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
import gl2lib as G


def In_all(N, H, lnn, h=0.05, chunk=200000):
    lo, hi = H.support()
    t = np.arange(lo, hi + h / 2, h)
    F = H(t) * np.exp(2j * G.theta(N, t)) * h
    F[0] *= 0.5; F[-1] *= 0.5
    out = np.empty(len(lnn), complex)
    for s in range(0, len(lnn), chunk):
        x = lnn[s:s + chunk]
        z = np.exp(-1j * h * x)
        acc = np.zeros(len(x), complex)
        for Fj in F[::-1]:                      # Horner: sum_j F_j z^j
            acc = acc * z + Fj
        out[s:s + chunk] = acc * np.exp(-1j * t[0] * x)
    return out


def run(cid, T1, T2, D, pari, out, do_voronoi=True):
    C = G.load_curves()[cid]; N, w = C['N'], C['w']
    gam = G.load_zeros(cid)
    H = G.Weight('box', T1, T2, D)
    lo, hi = H.support()
    Z = G.zero_side(N, gam, H)
    xlo = N * lo ** 2 / (4 * math.pi ** 2) * math.exp(-1.6); xhi = N * hi ** 2 / (4 * math.pi ** 2) * math.exp(1.6)
    t0 = time.time()
    n, Lam, _ = G.prime_power_data(pari, C['ainv'], C['bad'], int(xhi))
    m = n >= xlo
    I = In_all(N, H, np.log(n[m]))
    P = -(1 / G.TWO_PI) * np.sum(Lam[m] * n[m] ** -0.5 * I)
    # the same with half the step and a wider n-range (convergence)
    out('%s N=%d w=%+d box [%g,%g] Delta=%g: prime powers %d..%d (%d terms, %.0fs)' % (cid, N, w, T1, T2, D, int(xlo), int(xhi), m.sum(), time.time() - t0))
    out('   (A) zeros %.6f %+.6fi | primes (exact I_n) %.6f %+.6fi | |diff| %.2e' % (Z.real, Z.imag, P.real, P.imag, abs(Z - P)))
    Ps = G.prime_side(N, n, Lam, H, mode='exact')
    Ps2 = G.prime_side(N, n, Lam, H, mode='exact', second_order=True)
    out('       stationary phase: leading %.4f %+.4fi (|diff to exact| %.3f), with 2nd order %.4f %+.4fi (%.3f)' % (
        Ps.real, Ps.imag, abs(Ps - P), Ps2.real, Ps2.imag, abs(Ps2 - P)))
    if not do_voronoi:
        return
    t0 = time.time()
    nn, lam = G.lambda_all(pari, C['ainv'], int(xhi))
    mm = nn >= xlo
    I2 = In_all(N, H, np.log(nn[mm]))
    left = np.sum(lam[mm] * nn[mm] ** -0.5 * I2)
    # dual side: Hhat(-ln m) = int H(t) e^{i t ln m} dt, exact for the Gaussian-edged box:
    #   = (e^{i T2 u} - e^{i T1 u})/(i u) * exp(-Delta^2 u^2/4), u = ln m;  m = 1: T2 - T1
    mvals = np.arange(1, 50, dtype=float); u = np.log(mvals)
    Hh = np.empty(len(u), complex); Hh[0] = T2 - T1
    Hh[1:] = (np.exp(1j * T2 * u[1:]) - np.exp(1j * T1 * u[1:])) / (1j * u[1:]) * np.exp(-D ** 2 * u[1:] ** 2 / 4)
    right = w * np.sum(lam[:49] * mvals ** -0.5 * Hh)
    out('   (B) Voronoi: sum lambda(n) n^-1/2 I_n = %.6f %+.6fi | w sum_m lambda(m) m^-1/2 Hhat(-ln m) = %.6f %+.6fi (m=1 term %.3f) | |diff| %.2e  (%d n, %.0fs)' % (
        left.real, left.imag, right.real, right.imag, w * (T2 - T1), abs(left - right), mm.sum(), time.time() - t0))


if __name__ == '__main__':
    pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6)
    f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'exact_check_out.txt'), 'a')
    def out(s):
        print(s); f.write(s + '\n'); f.flush()
    cid, T1, T2, D = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    run(cid, T1, T2, D, pari, out, do_voronoi=(len(sys.argv) < 6 or sys.argv[5] != 'novor'))
