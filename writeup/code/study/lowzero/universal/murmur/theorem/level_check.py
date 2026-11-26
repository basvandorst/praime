"""Checks for level_theorem.md (fixed squarefree level N; own code; PARI via cypari2 for modular forms and class numbers).

Sections (labels as in the note):
 A  Euler-product identities: the dimension-weighted average over squarefree N of the per-level factors = Zubrilina's
    constants (exact, per prime); brute-force averages of phi(N) M^(N)(y) over squarefree N <= X against M_zub.
 B  local means E_l^(N)(r) over p in Z_l^x by exact enumeration: 1 if l | N (l = 2 included), generic otherwise.
 C  trace formula of T_p W_N on S_k^new(N), N squarefree composite, against PARI trace forms; old forms give trace 0;
    eps_f = (-1)^(k/2) w_N(f) = root number; dimension of the new space.
 D  the densities M^(N): (a) prime level, (c) the table for N = q, 2q, 3q, 6q, 5q.
 E  discovery levels near 1e4 (seen-type data, declared): window means vs M^(N), M_prime, M_zub, the other types;
    class-number means vs nu_N(r); calibration of the errors by the scatter BETWEEN levels; writes the frozen
    predictions (blind levels) and the dry-run predictions (discovery levels).
Usage: level_check.py [sections]   (default ABCDE)
"""
import sys, math, json, hashlib, os
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import cypari2

pari = cypari2.Pari()
pari.allocatemem(300 * 10 ** 6)
HERE = os.path.dirname(os.path.abspath(__file__))
NCHK = NFAIL = 0
OUT = {}


def check(name, ok, info=''):
    global NCHK, NFAIL
    NCHK += 1
    NFAIL += (not ok)
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)


# ------------------------------------------------------------------ densities
PR = [int(p) for p in sp.primerange(2, 10 ** 6)]


def u_(l): return 1.0 / ((l - 1) ** 2 * (l + 1))


C0 = math.exp(sum(math.log1p(-u_(l)) for l in PR))          # prod_l (1 - u_l); tail beyond 1e6 < 1e-12
Z2 = math.pi ** 2 / 6


def e_loc(l, r):
    """generic local mean E_l(r) (l not dividing N): l^2/(l^2-1) if l | r (r = 0 included), 1 - u_l otherwise."""
    return l * l / (l * l - 1.0) if r % l == 0 else 1.0 - u_(l)


def nu_prime(r):
    v = C0
    for l in sp.factorint(r):
        v *= (l * l / (l * l - 1.0)) / (1.0 - u_(l))
    return v


def nu_N(r, S):
    """nu_N(r) = prod_{l not| N} E_l(r) = nu'(r) / prod_{l in S} E_l(r), S = the primes dividing N."""
    v = nu_prime(r)
    for l in S:
        v /= e_loc(l, r)
    return v


def kappa(S): return float(np.prod([l / (l - 1.0) for l in S])) if S else 1.0
def E0(S): return Z2 * (float(np.prod([1 - 1.0 / l ** 2 for l in S])) if S else 1.0)


def Ucheb(n, t):
    a, b = 1.0, 2 * t
    if n == 0: return 1.0
    for _ in range(n - 1):
        a, b = b, 2 * t * b - a
    return b


def M_level(y, S, k=2):
    """M_k^(N)(y), N squarefree with prime set S (S = () : the generic density M_k of Theorem 2)."""
    sg = (-1) ** (k // 2 - 1)
    out = (12 / math.pi) * E0(S) * math.sqrt(y)
    r = 1
    while r * r < 4 * y:
        out += sg * (12 / math.pi) * nu_N(r, S) * Ucheb(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    out = kappa(S) * out / (k - 1)
    if k == 2:
        out -= 12 * kappa(S) * y
    return out


def _ep(fn): return math.exp(sum(math.log(fn(float(p))) for p in PR))


ZA = 2 * math.pi * _ep(lambda p: (1 - p - 2 * p ** 2 + p ** 4) / (p ** 4 - 2 * p ** 2 + p))
ZB = 2 * math.pi * _ep(lambda p: (-1 + p ** 2 + p ** 3) / (p * (-1 + p + p ** 2)))
ZG = 12 * _ep(lambda p: p * (1 + p) / (-1 + p + p ** 2))


def nu_zub(r):
    v = 1.0
    for l in sp.factorint(r):
        v *= 1 + l * l / (l ** 4 - 2 * l * l - l + 1)
    return v


def M_zub(y, k=2):
    """Zubrilina, arXiv 2310.07681 Thm 1 (with U_{k-2} inside the r-sum, as corrected in arXiv 2609.29631)."""
    out = ZB * math.sqrt(y) / (k - 1) - (ZG * y if k == 2 else 0.0)
    r = 1
    while r * r < 4 * y:
        out += ZA * (-1) ** (k // 2 - 1) / (k - 1) * nu_zub(r) * Ucheb(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    return out


# ------------------------------------------------------------------ A
def sec_A():
    bad = 0
    for l in sp.primerange(2, 1000):
        L = sp.Integer(l)
        e0 = L ** 2 / (L ** 2 - 1); e1 = 1 - 1 / ((L - 1) ** 2 * (L + 1))
        avg = lambda e: (L * e + 1) / (L + 1)      # P(l | N) = 1/(l+1) among squarefree N; factor 1 at l | N
        bad += sp.simplify(avg(e0) - (1 + L / ((L + 1) ** 2 * (L - 1)))) != 0                 # Zubrilina's A (r = 0)
        bad += sp.simplify(avg(e1) - (L ** 4 - 2 * L ** 2 - L + 1) / (L ** 2 - 1) ** 2) != 0    # Zubrilina's B (l !| r)
        bad += sp.simplify(avg(e0) / avg(e1) - (1 + L ** 2 / (L ** 4 - 2 * L ** 2 - L + 1))) != 0  # nu(r) factor (l | r)
        phi_avg = (L + 1 - 1 / L) / (L + 1)                                                    # mean of phi(N)/N, local
        e0_avg = (L + 1 - L ** -2) / (L + 1)                                                    # mean of prod_{q|N}(1-q^-2)
        bad += sp.simplify(e0_avg / phi_avg - (L ** 3 + L ** 2 - 1) / (L * (L ** 2 + L - 1))) != 0  # beta / 2 pi
        bad += sp.simplify(1 / phi_avg - L * (1 + L) / (L ** 2 + L - 1)) != 0                   # gamma / 12
        bad += sp.simplify((1 - L ** -2) * avg(e1) / phi_avg - (1 - L - 2 * L ** 2 + L ** 4) / (L ** 4 - 2 * L ** 2 + L)) != 0  # alpha
    check('A1 dimension-weighted average over squarefree N of the per-level local factors = Zubrilina\'s A, B, nu(r), '
          'alpha, beta, gamma (exact, every prime l < 1000)', bad == 0)
    # brute force: average of phi(N) M^(N)(y) over squarefree N <= X, divided by the average of phi(N)
    X = 300000
    Ns = np.arange(1, X + 1)
    sqf = np.ones(X + 1, bool)
    for l in sp.primerange(2, int(math.isqrt(X)) + 1):
        sqf[l * l::l * l] = False
    ys = (0.1, 0.5, 1.0, 2.0, 3.3)
    # multiplicative per-N factors via a sieve over primes: phi/N, prod(1-q^-2), and prod_{q|N} 1/E_q(r) for r = 1, 2, 3
    fphi = np.ones(X + 1); f0 = np.ones(X + 1); fr = {r: np.ones(X + 1) for r in (1, 2, 3)}
    for l in sp.primerange(2, X + 1):
        fphi[l::l] *= 1 - 1.0 / l; f0[l::l] *= 1 - 1.0 / l ** 2
        for r in (1, 2, 3):
            fr[r][l::l] /= e_loc(l, r)
    sel = sqf[1:]
    phiN = (Ns * fphi[1:])[sel]; Nn = Ns[sel].astype(float)
    worst = 0.0; rows = []
    for y in ys:
        # phi(N) M^(N)(y) = N [2 pi sqrt(y) prod(1-q^-2) + (12/pi) sum_r nu'(r) prod 1/E_q(r) sqrt(4y-r^2)] - 12 N y
        tot = Nn * (2 * math.pi * math.sqrt(y) * f0[1:][sel] - 12 * y)
        r = 1
        while r * r < 4 * y:
            tot = tot + Nn * (12 / math.pi) * nu_prime(r) * fr[r][1:][sel] * math.sqrt(4 * y - r * r)
            r += 1
        avg = tot.sum() / phiN.sum()
        worst = max(worst, abs(avg - M_zub(y)))
        rows.append('%.1f: %.4f vs %.4f' % (y, avg, M_zub(y)))
    check('A2 sum_{N<=3e5 sqfree} phi(N) M^(N)(y) / sum phi(N) = M_zub(y) (Zubrilina) at y = 0.1, 0.5, 1, 2, 3.3',
          worst < 5e-3, '; '.join(rows))


# ------------------------------------------------------------------ B
def local_F_bounds(Dl, l, J):
    """Delta = -D mod l^J for each sample; returns (lo, hi) enclosure of F_l(D) = sum_{j<a} l^-j + l^-a/(1 - chi(l)/l)."""
    x = 1.0 / l
    S = lambda a: (1 - x ** a) / (1 - x)
    v = np.zeros(len(Dl), dtype=np.int64)
    for kk in range(1, J + 1):
        v += (Dl % l ** kk == 0)
    vv = np.minimum(v, J)
    unit = Dl // (l ** vv)
    if l > 2:
        und = v >= J
        leg = np.array([0] + [int(sp.legendre_symbol(t, l)) for t in range(1, l)])
        a = vv // 2
        chi = np.where(vv % 2 == 1, 0, leg[unit % l])
        a0 = J // 2
    else:
        und = v >= J - 3
        uu = unit % 8
        odd = vv % 2 == 1
        a = np.where(odd, (vv - 3) // 2, np.where(uu % 4 == 1, vv // 2, vv // 2 - 1))
        chi = np.where(odd | (uu % 4 == 3), 0, np.where(uu == 1, 1, -1))
        a0 = max((J - 7) // 2, 0)
    a = np.maximum(a, 0)
    F = S(a) + x ** a / (1 - chi * x)
    lo = np.where(und, S(a0), F)
    hi = np.where(und, 1 / (1 - x), F)
    return lo, hi


def sec_B():
    J = {2: 19, 3: 11, 5: 7, 7: 6}
    levels = (10007, 10009, 10006, 9993, 10014, 10015, 10031, 15, 21, 35, 105, 210, 30030)
    nc = 0; okall = True; worst = 0.0; bad = []
    for l in (2, 3, 5, 7):
        q = l ** J[l]
        p = np.arange(q, dtype=np.int64); p = p[p % l != 0]
        for N in levels:
            for r in range(0, 13):
                Dl = (r * r * N * N - 4 * p * N) % q
                lo, hi = local_F_bounds(Dl, l, J[l])
                pred = 1.0 if N % l == 0 else e_loc(l, r)
                ok = lo.mean() - 1e-12 <= pred <= hi.mean() + 1e-12 and hi.mean() - lo.mean() < 2e-3
                okall &= ok; nc += 1; worst = max(worst, hi.mean() - lo.mean())
                if not ok: bad.append((l, N, r, lo.mean(), hi.mean(), pred))
    check('B1 local means E_l^(N)(r) over p in Z_l^x: = 1 if l | N (l = 2, 3, 5, 7), = generic E_l(r) if l !| N',
          okall, '%d cases (13 levels: odd N in all classes mod 8, N even, 3|N, 5|N, 7|N; r = 0..12), enclosure width <= %.1e %s'
          % (nc, worst, bad[:3]))


# ------------------------------------------------------------------ C
def H(D): return Fr(str(pari('qfbhclassno(%d)' % D)))


def U_t2(n, t2):
    """U_n(t) for even n as a polynomial in t^2 (exact)."""
    return {0: 1, 2: 4 * t2 - 1, 4: 16 * t2 * t2 - 12 * t2 + 1}[n]


def trace_rhs(p, N, k):
    s = H(4 * p * N) / 2
    r = 1
    while r * r * N < 4 * p:
        s += (-1) ** (k // 2 - 1) * U_t2(k - 2, Fr(r * r * N, 4 * p)) * H(4 * p * N - r * r * N * N)
        r += 1
    return s - (p + 1 if k == 2 else 0)


def omega(N): return len(sp.factorint(N))


def sec_C():
    tot = 0
    for k, levels, pmax in ((2, (6, 10, 14, 15, 21, 22, 26, 30, 33, 34, 35, 38, 39, 42, 66, 70, 78, 105), 300),
                            (4, (6, 10, 14, 15, 21, 30, 42), 160), (6, (6, 10, 15), 100)):
        bad = n = 0
        for N in levels:
            c = pari('mfcoefs(mftraceform([%d, %d], 0), %d)' % (N, k, N * pmax))
            sgn = (-1) ** (k // 2 + omega(N))
            for p in sp.primerange(2, pmax + 1):
                if N % p == 0:
                    continue
                lhs = sgn * Fr(str(c[N * p])) / Fr(N * p) ** (k // 2 - 1)   # = sum_f eps_f sqrt(p) lambda_f(p)
                bad += lhs != trace_rhs(p, N, k); n += 1
        tot += n
        check('C%d trace formula of T_p W_N on S_%d^new(N), squarefree N in %s, all primes p <= %d, p !| N (p = 2 included)'
              % (k // 2, k, levels, pmax), bad == 0, '%d pairs, exact; eps_f = (-1)^(k/2+omega(N)) a_f(N)/N^(k/2-1)' % n)
    # old forms: Tr(T_p W_N) on the full cusp space = on the new space
    bad = n = 0
    for N, k in ((30, 2), (42, 2), (66, 2), (70, 2), (30, 4), (42, 4)):
        trs = []
        for fl in (0, 1):
            pari('mf=mfinit([%d,%d],%d); A=mfatkininit(mf,%d)' % (N, k, fl, N))
            trs.append([pari('trace(mfheckemat(mf,%d)*A[2]/A[3])' % p) for p in sp.primerange(7, 60) if N % p])
        bad += trs[0] != trs[1]; n += len(trs[0])
    check('C4 old forms contribute 0: Tr(T_p W_N | S_k(N)) = Tr(T_p W_N | S_k^new(N)), N = 30, 42, 66, 70 (k=2), 30, 42 (k=4)',
          bad == 0, '%d traces' % n)
    # root numbers and Atkin-Lehner signs
    bad = n = 0
    for N, k in ((30, 2), (66, 2), (70, 2), (30, 4), (42, 4), (15, 6), (21, 4)):
        pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf)' % (N, k))
        for i in range(1, int(pari('#B')) + 1):
            if int(pari('poldegree(mffields(mf)[%d])' % i)) != 1:
                continue
            eps = int(pari('lfunrootres(lfunmf(mf,B[%d]))[3]' % i))
            wN = int(pari('mfatkineigenvalues(mf,%d)[%d][1]' % (N, i)))
            bad += eps != (-1) ** (k // 2) * wN
            for q in sp.factorint(N):
                wq = int(pari('mfatkineigenvalues(mf,%d)[%d][1]' % (q, i)))
                aq = int(pari('mfcoef(B[%d],%d)' % (i, q)))
                bad += aq != -wq * q ** (k // 2 - 1)
            n += 1
    check('C5 root number (lfunrootres) = (-1)^(k/2) w_N(f) and a_f(q) = -q^(k/2-1) w_q(f), rational newforms', bad == 0,
          '%d forms at 7 (level, weight) pairs' % n)
    worst = {2: 0.0, 4: 0.0}; bad = 0
    for N in range(2, 3001):
        if not sp.factorint(N) or max(sp.factorint(N).values()) > 1:
            continue
        ph = int(sp.totient(N)); w = omega(N)
        for k in (2, 4):
            dm = int(pari('mfdim([%d,%d],0)' % (N, k)))
            dev = abs(dm - Fr((k - 1) * ph, 12))
            worst[k] = max(worst[k], float(dev) / 2 ** w)
            bad += dev > 2 ** (w + 1)
    check('C6 |dim S_k^new(N) - (k-1) phi(N)/12| <= 2^(omega(N)+1), squarefree N <= 3000, k = 2, 4', bad == 0,
          'max |dev|/2^omega: k=2 %.3f, k=4 %.3f' % (worst[2], worst[4]))


# ------------------------------------------------------------------ D
TYPES = {'q': (), '2q': (2,), '3q': (3,), '6q': (2, 3), '5q': (5,)}


def sec_D():
    q = 1000003
    worst = max(abs(M_level(y, (q,), k) - M_level(y, (), k)) for k in (2, 4) for y in np.linspace(0.01, 6, 120))
    check('D1 (a) N = q prime: M^(q) = M_k (Theorem 2) + O(1/q): max dev at q = 1000003, k = 2, 4, y <= 6', worst < 1e-4,
          '%.1e' % worst)
    ys = (0.1, 0.2, 0.5, 1.0, 2.0)
    tab = {t: [M_level(y, S) for y in ys] for t, S in TYPES.items()}
    tab['zub'] = [M_zub(y) for y in ys]
    print('     (c) k = 2, q -> infinity:  y = ' + '  '.join('%7.2f' % y for y in ys))
    for t in tab:
        print('     %-4s ' % t + ' '.join('%8.4f' % v for v in tab[t]))
    tab4 = {t: [M_level(y, S, 4) for y in ys] for t, S in TYPES.items()}
    tab4['zub'] = [M_zub(y, 4) for y in ys]
    print('     (c) k = 4:')
    for t in tab4:
        print('     %-4s ' % t + ' '.join('%8.4f' % v for v in tab4[t]))
    OUT['table'] = {'y': ys, 'k2': tab, 'k4': tab4}
    # explicit constants
    for t, S in TYPES.items():
        print('     %-3s kappa = %.4f, prod(1+1/l) = %.4f, nu_N(1..4) = %s' % (t, kappa(S), kappa(S) * E0(S) / Z2,
              ', '.join('%.5f' % nu_N(r, S) for r in (1, 2, 3, 4))))
    check('D2 nu\'(1) = C0 = %.9f; E0 = zeta(2) prod(1-q^-2); kappa E0 12/pi = 2 pi prod(1+1/q)' % C0,
          all(abs(kappa(S) * E0(S) * 12 / math.pi - 2 * math.pi * np.prod([1 + 1.0 / l for l in S])) < 1e-12 for S in TYPES.values()))


# ------------------------------------------------------------------ E (discovery) and the frozen predictions
EDGES = np.arange(0, 4.0 + 1e-9, 0.25)
GP_LEVEL = ('lvl(N,Y)=my(v=List());forprime(p=2,Y*N,if(N%p,my(h=vector(4),r=0,s=0);while(r^2*N<4*p,'
            'my(c=qfbhclassno(4*p*N-r^2*N^2));if(r<4,h[r+1]=c);s+=if(r,c,c/2);r++);'
            'listput(v,[p,h[1],h[2],h[3],h[4],s-(p+1)])));Vec(v)')
pari(GP_LEVEL)


def level_data(N):
    """rows [p, H(4pN), H(4pN - N^2), H(4pN - 4N^2), H(4pN - 9N^2), T(p)] for primes p <= 4N, p !| N, with
    T(p) = H(4pN)/2 + sum_{r>=1, r^2 N < 4p} H(4pN - r^2 N^2) - (p + 1) = sum_f eps_f a_f(p) (k = 2, trace formula C1)."""
    m = pari('lvl(%d,4)' % N)
    return np.array([[float(x) for x in row] for row in m])


def model_means(N, S, ps, fn):
    y = ps / N; w = np.log(ps)
    vals = np.array([fn(t) for t in y])
    out = []
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        sel = (y > a) & (y <= b)
        out.append(float(np.sum(w[sel] * vals[sel]) / np.sum(w[sel])))
    return out


def own_exact(N, S, dim):
    """per-prime main term with the exact dimension and the exact (p+1): (1/dim)[(N/pi)(E0 sqrt y + sum nu_N sqrt(4y-r^2)) - (p+1)]."""
    def f(y):
        s = E0(S) * math.sqrt(y)
        r = 1
        while r * r < 4 * y:
            s += nu_N(r, S) * math.sqrt(4 * y - r * r); r += 1
        return (N / math.pi * s - (y * N + 1)) / dim
    return f


def levels_of(m, thr, n):
    out = []; q = thr
    while len(out) < n:
        q = int(sp.nextprime(q))
        if m % q:
            out.append(m * q)
    return out


MODELS = {'q': (), '2q': (2,), '3q': (3,), '6q': (2, 3)}
MULT = {'q': 1, '2q': 2, '3q': 3, '6q': 6}


def predictions(levels_by_type):
    pred = {}
    for t, Ns in levels_by_type.items():
        for N in Ns:
            S = tuple(sorted(sp.factorint(N)))
            dim = int(pari('mfdim([%d,2],0)' % N))
            ps = np.array([int(p) for p in pari('primes([2, %d])' % (4 * N)) if N % int(p)], dtype=float)
            ent = {'type': t, 'N': N, 'dim': dim, 'nprimes': int(len(ps)), 'kappa': kappa(S),
                   'own': model_means(N, S, ps, own_exact(N, S, dim)),
                   'zub': model_means(N, S, ps, M_zub)}
            for t2, S2 in MODELS.items():
                ent['type_' + t2] = model_means(N, S, ps, lambda y, S2=S2: M_level(y, S2))
            pred[str(N)] = ent
    return pred


def window_means(N, M):
    ps = M[:, 0]; y = ps / N; w = np.log(ps)
    dim = int(pari('mfdim([%d,2],0)' % N))
    E = M[:, 5] / dim
    out = []
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        sel = (y > a) & (y <= b)
        out.append(float(np.sum(w[sel] * E[sel]) / np.sum(w[sel])))
    return np.array(out)


DISC = {'q': levels_of(1, 10006, 5), '2q': levels_of(2, 5002, 5), '3q': levels_of(3, 3330, 5), '6q': levels_of(6, 1668, 5)}
BLIND = {'q': levels_of(1, 30030, 10), '2q': levels_of(2, 15015, 10), '3q': levels_of(3, 10010, 10), '6q': levels_of(6, 5005, 10)}
SEEN_NEAR_3E4 = (30011, 30029)          # prime levels computed by the referee of mprime_theorem.md


def sec_E():
    allN = [N for Ns in BLIND.values() for N in Ns]
    check('E0 blind levels: 40 distinct squarefree levels, none seen (referee: 30011, 30029; round 160: N <= 3200; '
          'discovery: N <= 10254)', len(set(allN)) == 40 and not set(allN) & set(SEEN_NEAR_3E4)
          and all(max(sp.factorint(N).values()) == 1 for N in allN), 'range %d..%d' % (min(allN), max(allN)))
    predD = predictions(DISC)
    data = {}; zs = {t: [] for t in DISC}
    for t, Ns in DISC.items():
        for N in Ns:
            M = level_data(N)
            data[N] = window_means(N, M)
            S = tuple(sorted(sp.factorint(N)))
            ps = M[:, 0]
            for r in range(0, 4):
                D = 4 * ps * N - r * r * N * N
                sel = D > 0
                Lh = math.pi * M[sel, 1 + r] / np.sqrt(D[sel])
                pr = E0(S) if r == 0 else nu_N(r, S)
                zs[t].append((r, Lh.mean(), pr, (Lh.mean() - pr) / (Lh.std() / math.sqrt(len(Lh)))))
    # pooled class-number means per type and r (direct test of the local factors at the primes dividing N)
    worstz = 0.0
    for t in DISC:
        line = []
        for r in range(4):
            v = [z for z in zs[t] if z[0] == r]
            m = np.mean([z[1] for z in v]); pr = v[0][2]
            zz = np.mean([z[3] for z in v]) * math.sqrt(len(v))
            worstz = max(worstz, abs(zz)); line.append('r=%d %.4f/%.4f (z %.2f)' % (r, m, pr, zz))
        print('     %-3s mean pi H(D_r)/sqrt(D_r) / prediction: ' % t + '; '.join(line))
    check('E1 class-number means over primes = E0(N), nu_N(1..3) (5 levels per type pooled; |z| < 4)', worstz < 4, 'max |z| %.2f' % worstz)
    # calibration: residuals R = data - own prediction; errors from the scatter BETWEEN levels (not single values):
    # per type a scale c_T (rms of R over its 5 levels x 16 windows), a shape s_w pooled over the 20 levels (R/c_T),
    # floored at a quarter of its median; sigma_{T,w} = c_T s_w
    Rr = {N: data[N] - np.array(predD[str(N)]['own']) for t in DISC for N in DISC[t]}
    cT = {t: float(np.sqrt(np.mean([Rr[N] ** 2 for N in DISC[t]]))) for t in DISC}
    s_w = np.sqrt(np.mean([(Rr[N] / cT[t]) ** 2 for t in DISC for N in DISC[t]], axis=0))
    s_reg = np.maximum(s_w, 0.25 * float(np.median(s_w)))
    SIG = {t: cT[t] * s_reg for t in DISC}
    for t in DISC:
        print('     %-3s c_T = %.4f (c_T/kappa = %.4f), mean residual %.4f' % (t, cT[t], cT[t] / kappa(MODELS[t]),
              float(np.mean([Rr[N] for N in DISC[t]]))))
    print('     shape s_w (20 levels) = ' + ' '.join('%.3f' % v for v in s_w))
    # single-value scatter for comparison (as in the earlier notes): standard error of a window mean from the per-prime spread
    Nx = DISC['2q'][0]; Mx = level_data(Nx); ps = Mx[:, 0]; y = ps / Nx; w = np.log(ps)
    dimx = int(pari('mfdim([%d,2],0)' % Nx)); Ex = Mx[:, 5] / dimx
    se = []
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        sel = (y > a) & (y <= b); x = Ex[sel]; ww = w[sel]; m = np.sum(ww * x) / np.sum(ww)
        se.append(math.sqrt(np.sum(ww ** 2 * (x - m) ** 2)) / np.sum(ww))
    print('     single-value s.e. at N = %d: ' % Nx + ' '.join('%.3f' % v for v in se))
    print('     between-level sigma (2q):    ' + ' '.join('%.3f' % v for v in SIG['2q']))
    print('     ratio between-level / single-value (median over windows): %.2f' % float(np.median(SIG['2q'] / np.array(se))))
    # discovery chi2 (in-sample: the errors come from the same levels)
    res = {}
    for t, Ns in DISC.items():
        res[t] = {}
        for mdl in ['own', 'zub'] + ['type_' + x for x in MODELS]:
            c2 = 0.0
            for N in Ns:
                sig = SIG[t]
                c2 += float(np.sum(((data[N] - np.array(predD[str(N)][mdl])) / sig) ** 2))
            res[t][mdl] = c2
        print('     discovery %-3s chi2 on %d terms: ' % (t, 16 * len(Ns)) + ', '.join('%s %.1f' % (m, v) for m, v in res[t].items()))
    ok = all(res[t]['own'] < 2 * 16 * len(DISC[t]) and min(v - res[t]['own'] for m, v in res[t].items() if m not in ('own', 'type_' + t)) > 100 for t in DISC)
    check('E2 discovery (in-sample errors): each type fits its own M^(N) (chi2 < 2 x terms); every other density: delta chi2 > 100', ok)
    rmsd = {t: float(np.sqrt(np.mean([(data[N] - np.array(predD[str(N)]['own'])) ** 2 for N in DISC[t]]))) for t in DISC}
    print('     discovery rms(data - own) per type: ' + ', '.join('%s %.4f' % kv for kv in rmsd.items()))
    for t in DISC:
        N = DISC[t][0]
        print('     N = %d (%s): data ' % (N, t) + ' '.join('%.3f' % v for v in data[N]))
        print('     %s own   ' % (' ' * len(str(N))) + ' '.join('%.3f' % v for v in predD[str(N)]['own']))
    # freeze
    meta = {'windows': [[float(a), float(b)] for a, b in zip(EDGES[:-1], EDGES[1:])], 'weight': 2, 'ymax': 4,
            'weights': 'log p', 'sigma_rule': 'sigma_{T,w} = c_T s_w: c_T = rms of (data - own) over the 5 discovery levels of type T, '
            's_w = pooled between-level shape over the 20 discovery levels, floored at 1/4 of its median; not scaled with N',
            'sigma': {t: [float(v) for v in SIG[t]] for t in SIG},
            'decision': {'D1_own_chi2_max': 256.0, 'D2_alt_delta_chi2_min': 100.0, 'D2_alt_chi2_min': 256.0, 'terms_per_type': 160}}
    predB = predictions(BLIND)
    for nm, pr, lv in (('predictions_blind.json', predB, BLIND), ('predictions_dry.json', predD, DISC)):
        doc = {'meta': meta, 'levels': {t: lv[t] for t in lv}, 'pred': pr}
        txt = json.dumps(doc, indent=1, sort_keys=True)
        fn = os.path.join(HERE, nm)
        if os.path.exists(fn):          # frozen: never overwrite, only check that it is reproduced byte for byte
            old = open(fn).read()
            check('E3 %s reproduced byte for byte (sha256 %s)' % (nm, hashlib.sha256(old.encode()).hexdigest()), old == txt)
        else:
            with open(fn, 'w') as fh:
                fh.write(txt)
            print('     wrote %s sha256 %s' % (nm, hashlib.sha256(txt.encode()).hexdigest()))
    # expected separations on the blind levels (no data): sum (own - alt)^2 / sigma^2 per type
    for t, Ns in BLIND.items():
        line = []
        for mdl in ['zub'] + ['type_' + x for x in MODELS if x != t]:
            e = sum(float(np.sum(((np.array(predB[str(N)]['own']) - np.array(predB[str(N)][mdl])) / SIG[t]) ** 2)) for N in Ns)
            line.append('%s %.0f' % (mdl, e))
        print('     blind %-3s expected delta chi2 (sum of squared separations / sigma^2): ' % t + ', '.join(line))
    with open(os.path.join(HERE, 'discovery_data.json'), 'w') as fh:
        json.dump({str(N): [float(v) for v in data[N]] for N in data}, fh)
    sf = os.path.join(HERE, 'state_dry.json')
    if os.path.exists(sf):              # the dry run of blind_test.py (separate counting code) on the same 20 levels
        st = json.load(open(sf))
        dev = max(abs(x - y) for N in data for x, y in zip(st[str(N)]['mean'], data[N]))
        check('E4 blind_test.py dry run (own counting code, exact rationals) = this script on the 20 discovery levels',
              dev < 1e-9 and len(st) == len(data), 'max dev %.1e' % dev)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'ABCDE'
    for ch in which:
        globals()['sec_' + ch]()
    print('%d checks, %d failed' % (NCHK, NFAIL))
