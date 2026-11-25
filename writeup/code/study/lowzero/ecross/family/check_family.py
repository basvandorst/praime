"""check_family.py: checks every identity and number used in family_proof.md.

Run:  python check_family.py a   (Kloosterman/Ramanujan identities, Bessel-Mellin, Sato-Tate local algebra)
      python check_family.py b   (Petersson formula and its root-number form; harmonic weights; root numbers)
      python check_family.py c   (harmonic and natural averages of c_f over full Hecke bases at prime levels)
      python check_family.py d   (Delta: c_Delta, 1/omega_Delta, zeros to 420, one-point sums, exact explicit formula,
                                  prime side, Petersson c-decomposition of the resonance, resonance sums)
Each section prints its numbers and asserts the identity (PASS/FAIL line).
"""
import sys, math, time
sys.path.insert(0, '.')
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from famlib import (theta, kloosterman, kloosterman_vec_m, ramanujan, mobius, petersson_delta, rootnumber_petersson,
                    primes_upto, family_data, c_f, units_and_inverses, _phi, stat_weights)
from scipy.special import jv

RES = []


def check(name, ok, info=''):
    RES.append((name, ok))
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)


# =========================================================================================== a
def sec_a():
    rng = np.random.default_rng(1)
    # a1: S(Nn,1;Nm) = mu(N) S(n, Nbar; m) for N prime, (m,N)=1
    worst = 0.0
    for N in [5, 7, 11, 13, 37, 101]:
        for _ in range(12):
            m = int(rng.integers(1, 60))
            if m % N == 0:
                continue
            n = int(rng.integers(1, 10 ** 6))
            Nbar = pow(N, -1, m) if m > 1 else 0
            lhs = kloosterman(N * n, 1, N * m)
            rhs = -kloosterman(n, Nbar, m)
            worst = max(worst, abs(lhs - rhs))
    check('a1 S(Nn,1;Nm) = mu(N) S(n,Nbar;m)', worst < 1e-9, 'max err %.1e' % worst)
    # a2: sum_{a in (Z/m)*} S(l a, b; m) = mu(m) c_m(l)  for (b,m)=1
    worst = 0.0
    for m in range(1, 41):
        for l in [1, 2, 3, 4, 6, 9, 12, 25, 36]:
            for b in [1, 5, 7]:
                if math.gcd(b, m) != 1:
                    continue
                if m == 1:
                    s = 1.0
                else:
                    x, _ = units_and_inverses(m)
                    s = sum(kloosterman(l * int(a), b, m) for a in x)
                worst = max(worst, abs(s - mobius(m) * ramanujan(m, l)))
    check('a2 sum_a* S(la,b;m) = mu(m) c_m(l)', worst < 1e-8, 'max err %.1e' % worst)
    # a3: Parseval: (1/phi(m)) sum_a* |S(Na,1;Nm)|^2 <= m, hence sum_chi |c_chi| <= m
    worst = 0.0
    for N in [7, 37]:
        for m in range(2, 80):
            if m % N == 0:
                continue
            x, _ = units_and_inverses(m)
            v = kloosterman_vec_m(N * x, 1, N * m)
            worst = max(worst, np.mean(v ** 2) / m)
    check('a3 (1/phi(m)) sum |K(a)|^2 <= m', worst <= 1 + 1e-9, 'max ratio %.3f' % worst)
    # a4: Bessel-Mellin  int_0^inf y^{s-1} J_nu(y) dy = 2^{s-1} Gamma((nu+s)/2)/Gamma((nu-s)/2+1)
    import mpmath as mp
    mp.mp.dps = 20
    worst = 0.0
    for nu in [1, 3, 11]:
        for s in [mp.mpc(0, -2 * 7.3), mp.mpc(1, -2 * 4.1), mp.mpc(0.6, -2 * 2.5)]:
            f = lambda y: y ** (s - 1) * mp.besselj(nu, y)
            val = mp.quad(f, mp.linspace(0, 200, 401)) + mp.quadosc(f, [200, mp.inf], omega=1)
            ex = 2 ** (s - 1) * mp.gamma((nu + s) / 2) / mp.gamma((nu - s) / 2 + 1)
            worst = max(worst, float(abs(val - ex) / abs(ex)))
    check('a4 Bessel Mellin transform (s = -2it, 1-2it, 0.6-2it)', worst < 1e-8, 'max rel err %.1e' % worst)
    # a5: the u-form used in the note: int_0^inf u^{-1/2-it} J_{k-1}(4 pi sqrt(u/(N m^2))) du = (sqrt N m/2pi) e^{-2i theta_{Nm^2}(t)}
    worst = 0.0
    for (k, N, m, t) in [(2, 37, 1, 9.0), (4, 13, 2, 6.5), (12, 1, 1, 15.0), (2, 11, 3, 4.0)]:
        a = 4 * mp.pi / mp.sqrt(N * m * m)
        # substitute y = a sqrt u: u = y^2/a^2, du = 2y/a^2 dy
        f = lambda y: (y ** 2 / a ** 2) ** mp.mpc(-0.5, -t) * mp.besselj(k - 1, y) * 2 * y / a ** 2
        val = mp.quad(f, mp.linspace(0, 200, 401)) + mp.quadosc(f, [200, mp.inf], omega=1)
        th = t / 2 * mp.log(N * m * m) - t * mp.log(2 * mp.pi) + mp.im(mp.loggamma(k / 2 + 1j * t))
        ex = mp.sqrt(N) * m / (2 * mp.pi) * mp.exp(-2j * th)
        worst = max(worst, float(abs(val - ex) / abs(ex)))
    check('a5 int u^{-1/2-it} J_{k-1}(4pi sqrt(u/(Nm^2))) du = (sqrt(N) m/2pi) e^{-2i theta_{Nm^2}(t)}', worst < 1e-8,
          'max rel err %.1e' % worst)
    # a6: on Re w = 0 the Mellin kernel has modulus exactly 1/|a + i y|, a = (k-1)/2
    worst = 0.0
    for k in [2, 4, 12]:
        a = (k - 1) / 2
        for y in [-30.0, -3.0, 0.4, 7.0, 55.0]:
            r = abs(mp.gamma(a + 1j * y) / mp.gamma(a + 1 - 1j * y))
            worst = max(worst, float(abs(r * abs(mp.mpc(a, y)) - 1)))
    check('a6 |Gamma(a+iy)/Gamma(a+1-iy)| = 1/|a+iy|', worst < 1e-12, 'max err %.1e' % worst)
    # a7: Sato-Tate local algebra: the conjectured twisted line amplitudes average to the family amplitudes
    import sympy as sp
    th_ = sp.symbols('th', real=True)
    X = sp.symbols('X', positive=True)
    def U(j):
        return sp.chebyshevu(j, sp.cos(th_)) if j >= 0 else 0
    def EST(expr):
        return sp.simplify(sp.integrate(sp.expand(expr * 2 / sp.pi * sp.sin(th_) ** 2), (th_, 0, sp.pi)))
    p = sp.symbols('p', positive=True)
    lamp = U(1)
    mu = {0: 1, 1: -lamp, 2: 1}
    ok = True
    rows = []
    for e in range(0, 4):
        for v in range(-4, 5):
            S = 0
            for h in range(0, 3):
                mm = v + h
                if mm < 0:
                    continue
                S += U(mm) * mu[h] * p ** (-sp.Rational(mm + h, 2))
            lhs = sp.simplify(EST(U(e) * S) / (1 - 1 / p))
            # family: r = l/m^2 locally: v = e - 2j, j = 0 or 1 (mu(p^j) c_{p^j}(p^e)/(phi(p^j) p^{e/2}))
            if (e - v) % 2 == 0 and e - v >= 0:
                j = (e - v) // 2
                if j == 0:
                    rhs = p ** (-sp.Rational(e, 2))
                elif j == 1:
                    rhs = -((p - 1) if e >= 1 else -1) / (p - 1) * p ** (-sp.Rational(e, 2))
                else:
                    rhs = 0
            else:
                rhs = 0
            d = sp.simplify(lhs - rhs)
            rows.append((e, v, lhs))
            ok = ok and d == 0
    check('a7 E_ST[lambda(p^e) S_p(v)]/(1-1/p) = family local factor mu(p^j)c_{p^j}(p^e)/(phi(p^j)p^{e/2}), e<=3, |v|<=4',
          ok)
    # a8: harmonic and natural Sato-Tate means of the local factor of c_f
    m1 = EST(1 - U(2) / p)
    m2 = EST((1 - U(2) / p) * sum(U(2 * j) * p ** (-j) for j in range(0, 8)))
    m3 = EST((1 - U(2) / p) ** 2)
    check('a8 E_ST[1-lambda(p^2)/p]=1, E_ST[(1-lambda(p^2)/p) Z_p]=1-p^-2 (+O(p^-8)), E_ST[(1-lambda(p^2)/p)^2]=1+p^-2',
          sp.simplify(m1 - 1) == 0 and sp.simplify(m2 - (1 - p ** -2)) == 0
          and sp.simplify(m3 - 1 - p ** -2) == 0, 'means: %s | %s | %s' % (m1, sp.simplify(m2), m3))


# =========================================================================================== b
def sec_b():
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(800 * 10 ** 6)
    pari.set_real_precision(30)
    # b1: Delta (level 1, weight 12): omega = Delta(1,1) = Gamma(11)/((4pi)^11 <D,D>) = 2pi^2/(11 L(1,Sym^2 D))
    tau = [int(x) for x in pari('mfcoefs(mfDelta(),250)')]
    om = petersson_delta(12, 1, 1, 1, 400)
    fs = family_data(pari, 1, 12, 250)
    LS1 = fs[0]['LS1']
    nD = 1.035362056804320922e-6     # <Delta,Delta> over SL2(Z)\H (literature value, used only as a cross-check)
    om2 = math.gamma(11) / ((4 * math.pi) ** 11 * nD)
    print('  omega_Delta: Petersson series %.15f, from <D,D> %.15f, 2pi^2/(11 L(1,Sym^2)) %.15f, L(1,Sym^2 Delta) = %.15f'
          % (om, om2, fs[0]['omega'], LS1))
    check('b1 omega_Delta three ways', abs(om - om2) < 1e-10 and abs(om - fs[0]['omega']) < 1e-10)
    worst = 0.0
    for n in [2, 3, 4, 5, 6, 7, 10, 11, 13, 97, 101, 199, 211, 241]:
        d = petersson_delta(12, 1, 1, n, 4000)
        worst = max(worst, abs(d / om - tau[n] / n ** 5.5))
    check('b2 Petersson for Delta: omega lambda(n) = Delta_{12,1}(1,n)', worst < 1e-12, 'max err %.1e' % worst)
    # b3, b4: weight 4 at levels 7, 13; weight 2 at level 37 (two rational forms, root numbers -1 and +1)
    for (N, k, Mm, tol) in [(7, 4, 600, 1e-6), (13, 4, 600, 1e-6), (37, 2, 1500, 2e-2), (11, 2, 1500, 2e-2)]:
        t0 = time.time()
        fs = family_data(pari, N, k, 60)
        w1 = w2 = 0.0
        for n in ([1, 2, 3, 4, 5, 6, 10] if k > 2 else [1, 2, 3, 5]):
            lhs = sum(f['omega'] * f['lam'][n] for f in fs)
            rhs = petersson_delta(k, N, 1, n, N * Mm)
            lhs2 = sum(f['omega'] * f['eps'] * f['lam'][n] for f in fs)
            rhs2 = rootnumber_petersson(k, N, n, Mm)
            w1 = max(w1, abs(lhs - rhs))
            w2 = max(w2, abs(lhs2 - rhs2))
        print('  N=%d k=%d: forms (eps, L(1,Sym^2), omega) %s; max |Petersson err| %.1e, |root-number Petersson err| %.1e '
              '(m <= %d), %.0fs' % (N, k, [(f['eps'], round(f['LS1'], 6), round(f['omega'], 6)) for f in fs], w1, w2, Mm,
                                    time.time() - t0))
        check('b3 N=%d k=%d: sum_f omega_f lam(n) = Delta(1,n) with omega_f = 2pi^2/((k-1) N L(1,Sym^2 f))' % (N, k), w1 < tol)
        check('b4 N=%d k=%d: sum_f omega_f eps_f lam(n) = -2pi N^-1/2 sum_m S(Nn,1;Nm)/m J_{k-1}(4pi sqrt(n/N)/m)' % (N, k),
              w2 < tol)
    # b5: eps_f = i^k mu(N) lambda_f(N) sqrt N against PARI's root number of L(s,f)
    ok = True
    for (N, k) in [(13, 4), (37, 2), (11, 4), (23, 2)]:
        pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
        nb = int(pari('#B'))
        for j in range(1, nb + 1):
            Ls = pari('L=lfunmf(mf,B[%d]); if(type(L[1])=="t_VEC" && #L>1 && type(L[1][1])!="t_INT", L, [L])' % j)
            aN = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, N))
            vals = [aN] if str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) != 't_VEC' else list(aN)
            for idx, e in enumerate(vals):
                lamN = float(complex(e[N]).real) / N ** ((k - 1) / 2)
                epsf = (-1) ** (k // 2) * (-1) * lamN * math.sqrt(N)
                try:
                    r = complex(pari('lfunrootres(%s)[3]' % ('L[%d]' % (idx + 1) if len(vals) > 1 else 'L')))
                except Exception:
                    r = complex(pari('lfunrootres(L)[3]'))
                ok = ok and abs(r - epsf) < 1e-6
                print('  N=%d k=%d form %d emb %d: i^k mu(N) lam(N) sqrt N = %+.6f, PARI root number %s' % (N, k, j, idx, epsf, r))
    check('b5 root number formula', ok)


# =========================================================================================== c
def sec_c():
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1000 * 10 ** 6)
    pari.set_real_precision(30)
    P = 5000
    rows = []
    for (k, Ns) in [(4, [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]),
                    (2, [11, 37, 43, 53, 61, 67, 73, 79, 83, 89, 97])]:
        for N in Ns:
            t0 = time.time()
            fs = family_data(pari, N, k, P)
            om = np.array([f['omega'] for f in fs])
            c = np.array([c_f(f['lam'], N, P) / f['LS1'] for f in fs])
            eps = np.array([f['eps'] for f in fs])
            sum_om = om.sum()
            Mh = float(np.sum(om * c))
            nat = float(np.mean(c))
            pred_nat = (6 / math.pi ** 2) / (1 - 1 / N)
            # orthogonality of eps_f and c_f (harmonic): sum omega eps c
            rows.append((k, N, len(fs), sum_om, Mh, 1 + 1 / N, nat, pred_nat, float(np.sum(om * eps * c)), float(np.sum(om * c * c))))
            print('  k=%d N=%3d #f=%2d  sum omega=%.4f  M_pred=sum^h c_f=%.4f (1+1/N=%.4f)  natural mean c_f=%.4f '
                  '(6/pi^2/(1-1/N)=%.4f)  sum^h eps c=%+.4f  sum^h c^2=%.4f  [%.1fs]' % (k, N, len(fs), sum_om, Mh, 1 + 1 / N,
                                                                                    nat, pred_nat, rows[-1][8], rows[-1][9],
                                                                                    time.time() - t0), flush=True)
    r = np.array(rows, float)
    sel = ((r[:, 0] == 4) & (r[:, 1] >= 13)) | ((r[:, 0] == 2) & (r[:, 1] >= 43))
    for kk in (4, 2):
        s_ = sel & (r[:, 0] == kk)
        print('  k=%d, %d levels: mean of sum^h c_f = %.3f (sd %.3f); mean of natural mean c_f = %.3f (limit 6/pi^2 = %.3f); '
              'mean sum omega = %.3f' % (kk, s_.sum(), r[s_, 4].mean(), r[s_, 4].std(), r[s_, 6].mean(), 6 / math.pi ** 2,
                                         r[s_, 3].mean()))
    check('c1 sum_f omega_f = Delta(1,1) close to 1 (|.-1| < 0.15) at every level N >= 13', np.all(np.abs(r[r[:, 1] >= 13, 3] - 1) < 0.15))
    check('c2 harmonic mean of c_f, averaged over the levels, within 0.15 of the limit 1 (small levels: computed, not a test '
          'of the limit)', abs(r[sel, 4].mean() - 1) < 0.15, 'mean %.3f' % r[sel, 4].mean())
    np.save('family_c_rows.npy', r)


# =========================================================================================== d
def sec_d():
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1000 * 10 ** 6)
    pari.set_real_precision(25)
    k, N = 12, 1
    P = 200000
    co = pari('mfcoefs(mfDelta(),%d)' % P)
    tau = np.array([float(int(co[i])) for i in range(P + 1)])
    nn = np.arange(P + 1, dtype=float)
    lam = np.zeros(P + 1)
    lam[1:] = tau[1:] / nn[1:] ** 5.5
    pr = primes_upto(P)
    e1 = lam[pr] ** 2 - 1
    p = pr.astype(float)
    B = math.exp(np.sum(np.log1p(-e1 / p) - np.log1p(-e1 / p + e1 / p ** 2 - 1 / p ** 3)))
    fs = family_data(pari, 1, 12, 10)
    LS1 = fs[0]['LS1']
    om = fs[0]['omega']
    cD = B / LS1
    direct = math.exp(np.sum(np.log1p(-e1 / p)))
    print('  Delta: L(1,Sym^2)=%.12f  B=%.10f  c_Delta=B/L=%.6f  direct prod_{p<=2e5}(1-lambda(p^2)/p)=%.6f  '
          '1/omega_Delta=%.6f' % (LS1, B, cD, direct, 1 / om))
    check('d1 c_Delta = B/L(1,Sym^2) agrees with the direct partial product (0.1%)', abs(direct / cD - 1) < 1e-3)
    t0 = time.time()
    g = np.array([float(x) for x in pari('LD=lfuncreate(lfunetaquo(Mat([1,24]))); lfunzeros(LD,420)')])
    print('  %d zeros of L(Delta) in (0,420] (%.1fs); counts vs theta/pi:' % (len(g), time.time() - t0),
          [(T, int(np.sum(g <= T)), round(float(theta(k, N, T)) / math.pi, 2)) for T in (100, 200, 300, 400)])
    th = lambda t: theta(k, N, t)
    out = []
    for (T1, T2) in [(20, 400), (20, 140), (140, 260), (260, 400), (60, 230), (230, 380)]:
        m = (g > T1) & (g <= T2)
        gg = g[m]
        sharp = np.sum(np.exp(2j * th(gg)))
        w = np.sin(math.pi * (gg - T1) / (T2 - T1)) ** 2
        hann = np.sum(w * np.exp(2j * th(gg)))
        law = -cD * (T2 - T1) / (2 * math.pi)
        out.append((T1, T2, sharp, hann, law))
        print('  window [%d,%d]: sharp %s  Hann %s | law -c T/2pi: %.2f (Hann %.2f) | ratios %.3f, %.3f | c=1 term only '
              '(1/omega): %.2f' % (T1, T2, np.round(sharp, 2), np.round(hann, 2), law, law / 2, sharp.real / law,
                                    hann.real / (law / 2), -(T2 - T1) / (2 * math.pi * om)))
    s0 = out[0]
    check('d2 Delta one-point sum on [20,400] within 5% of -c_Delta T/2pi (sharp and Hann)',
          abs(s0[2].real / s0[4] - 1) < 0.05 and abs(s0[3].real / (s0[4] / 2) - 1) < 0.05)
    check('d3 Delta one-point sum excludes the c=1 Petersson amplitude 1/omega (off by > 20 x noise)',
          abs(s0[2].real - (-(380) / (2 * math.pi * om))) > 100)
    # d4: exact explicit formula (Proposition 1, weight 12) on a Gaussian-smoothed box: zeros vs primes (completeness)
    from scipy.special import erf
    for (T1, T2, D) in [(100.0, 200.0, 5.0), (250.0, 330.0, 6.0)]:
        H = lambda t: 0.5 * (erf((t - T1) / D) - erf((t - T2) / D))
        zs = np.sum(H(g) * np.exp(2j * th(g)))
        tt = np.arange(T1 - 7 * D, T2 + 7 * D, 0.01)
        Ht = H(tt)
        ph = 2 * th(tt)
        thp = np.gradient(ph, tt)
        nlo = int(math.exp(thp.min() - 2.5))
        nhi = int(math.exp(thp.max() + 2.5)) + 1
        tot = 0.0 + 0.0j
        for pp in pr[pr <= nhi]:
            pp = int(pp)
            al = lam[pp]
            Uj = [1.0, al]
            j = 1
            q = pp
            while q <= nhi:
                if q >= nlo:
                    Lam = (Uj[j] - (Uj[j - 2] if j >= 2 else 0.0)) * math.log(pp)
                    In = np.trapezoid(Ht * np.exp(1j * (ph - tt * math.log(q))), tt)
                    tot += Lam * q ** -0.5 * In
                j += 1
                Uj.append(al * Uj[-1] - Uj[-2])
                q *= pp
        ps = -tot / (2 * math.pi)
        print('  exact explicit formula, box [%g,%g] Delta=%g: zeros %s  primes %s  |diff| %.2e' % (T1, T2, D, np.round(zs, 6),
                                                                                                  np.round(ps, 6), abs(zs - ps)))
        check('d4 Proposition 1 (weight 12) zeros = primes on box [%g,%g]' % (T1, T2), abs(zs - ps) < 1e-4)
    # d5: prime side of the one-point sum (leading stationary phase) on Hann [20,400]
    T1, T2 = 20.0, 400.0
    Hh = lambda t: np.where((t > T1) & (t < T2), np.sin(math.pi * (t - T1) / (T2 - T1)) ** 2, 0.0)
    xmax = int((T2 / (2 * math.pi)) ** 2 * 1.3) + 2
    prs = primes_upto(xmax)
    tp = 2 * math.pi * np.sqrt(prs)
    wp = stat_weights(k, N, prs, Hh)                  # exact stationary point, exact phase
    Ppr = np.sum(lam[prs] * np.log(prs) * wp)
    p2 = prs[prs ** 2 <= xmax]
    w2 = stat_weights(k, N, p2.astype(float) ** 2, Hh)
    Ppr2 = np.sum((lam[p2] ** 2 - 2) * np.log(p2) * w2)
    Zp = Ppr + Ppr2
    pref = -np.exp(1j * (2 * k - 1) * math.pi / 4) / math.sqrt(2)
    Zlead = pref * np.sum(lam[prs] * np.log(prs) * prs ** -0.25 * np.exp(-4j * math.pi * np.sqrt(prs)) * Hh(tp))
    print('  (simplified kernel t_n* = 2 pi sqrt n, phase -4 pi sqrt n + (k-1) pi/2: %s; its O(k^2/t) phase error is large '
          'for k = 12 at t < 60)' % np.round(Zlead, 3))
    print('  Hann [20,400]: zero side %s, prime side (leading stationary phase) %s, law %.2f' % (np.round(out[0][3], 3),
                                                                                                 np.round(Zp, 3), out[0][4] / 2))
    check('d5 Delta: zeros = leading-order prime side on Hann [20,400] within 0.5', abs(Zp - out[0][3]) < 0.5)
    # d6: Petersson c-decomposition of the prime side: lambda(p) = (2pi/omega) sum_c S(p,1;c)/c J_11(4pi sqrt p/c)
    cmax = 800
    contrib = []
    lam_rec = np.zeros(len(prs))
    for c in range(1, cmax + 1):
        if c == 1:
            Sc = np.ones(len(prs))
        else:
            Sa = kloosterman_vec_m(np.arange(c), 1, c)
            Sc = Sa[prs % c]
        term = 2 * math.pi / om * Sc / c * jv(11, 4 * math.pi * np.sqrt(prs) / c)
        lam_rec += term
        contrib.append(np.sum(term * np.log(prs) * wp))
    contrib = np.array(contrib)
    errrec = np.max(np.abs(lam_rec - lam[prs]))
    cum = np.cumsum(contrib)
    law_h = out[0][4] / 2
    print('  lambda(p) rebuilt from the Kloosterman series (c <= %d): max err %.1e' % (cmax, errrec))
    for C in [1, 2, 3, 5, 10, 15, 20, 30, 40, 60, 80, 120, 200, 800]:
        print('    sum_{c<=%3d}: %s   (as fraction of the law %.2f: %.3f)' % (C, np.round(cum[C - 1], 3), law_h,
                                                                               cum[C - 1].real / law_h))
    check('d6 Kloosterman series reproduces lambda_Delta(p), p <= %d' % xmax, errrec < 1e-9)
    check('d6b c = 1 term gives 1/omega of the amplitude (within 10%), all c give the zero side (within 0.5)',
          abs(contrib[0].real / (-(T2 - T1) / (4 * math.pi * om)) - 1) < 0.1 and abs(cum[-1] + Ppr2 - out[0][3]) < 0.5,
          'c=1: %.3f vs %.3f' % (contrib[0].real, -(T2 - T1) / (4 * math.pi * om)))
    # d7: resonance sums over primes and over integers (sharp), x up to 2e5
    for x in [5e4, 1e5, 2e5]:
        n = np.arange(1, int(x) + 1)
        Sint = np.sum(lam[1:int(x) + 1] * np.exp(-4j * math.pi * np.sqrt(n)))
        pp = pr[pr <= x]
        Spr = np.sum(lam[pp] * np.log(pp) * np.exp(-4j * math.pi * np.sqrt(pp)))
        q2 = pr[pr ** 2 <= x]
        Spr += np.sum((lam[q2] ** 2 - 2) * np.log(q2) * np.exp(-4j * math.pi * q2))
        A = 2 * math.sqrt(2) / 3 * np.exp(-1j * (2 * k - 1) * math.pi / 4) * x ** 0.75
        print('  x=%.0e: integers S/A = %s (theorem: 1);  prime powers S/A = %s (law: c_Delta = %.3f; c=1 term: %.3f)'
              % (x, np.round(Sint / A, 3), np.round(Spr / A, 3), cD, 1 / om))
    check('d7 integer resonance constant (2 sqrt2/3) e^{-i(2k-1)pi/4} x^{3/4} for k = 12 (|S/A-1| < 0.04, decreasing)', abs(Sint / A - 1) < 0.04)
    np.save('delta_out.npy', np.array([[o[0], o[1], o[2].real, o[2].imag, o[3].real, o[3].imag, o[4]] for o in out]))


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'a'
    t0 = time.time()
    {'a': sec_a, 'b': sec_b, 'c': sec_c, 'd': sec_d}[which]()
    print('section %s: %d/%d PASS, %.0fs' % (which, sum(1 for _, o in RES if o), len(RES), time.time() - t0))
