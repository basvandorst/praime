"""spectrometer.py: recover the conductor N, the root number and the constant c of an 'unknown' elliptic curve from its
list a(p), p <= X only (no zeros, no conductor given), with the chirp lines of spectro.md section 2.

1. lambda(p) = a(p)/sqrt p; spectrum S(alpha) = sum_p lambda(p) ln p Hann(sqrt p) e(-alpha sqrt p) on a fine grid
   (oversampling 32, alpha <= 3), p in [1e4, X].
2. For every candidate N' in [2, NMAX] (bad primes = primes dividing N', assumed multiplicative if a(l) = +-1, additive if 0):
   predicted relative line weights w_r = (c(r)/c) r^{-1/4} for r = m, m/2, m/3, m/4 (alpha_r = 2 sqrt(r/N') <= 3) from the
   a(l) of the data (good-prime local factors of lines_lib; bad multiplicative l: A(k) = lambda(l^k)(1 + 1/l), A(-1) =
   -lambda(l)/(l - 1), A(-2) = 0 (measured on integers for 11a, 37a: val_lines.py); additive: A(0) = l/(l-1), else 0).
   Matched filter: z(N', eps) = Re[eps e^{3 pi i/4} sum_r w_r S(alpha_r)] / (sigma sqrt(sum w_r^2)), sigma = sqrt(sum w^2/2).
3. Estimates: N_hat, eps_hat = argmax z; c_hat = slope / (unit amplitude), compared with the truth from PARI (ellglobalred,
   ellrootno, c = B/L(1, Sym^2) as in curves_c.py).
usage: spectrometer.py X NMAX"""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spectro_lib as S

X = int(float(sys.argv[1])); NMAX = int(sys.argv[2]); ONLY = sys.argv[3].split(',') if len(sys.argv) > 3 else None
P = S.pari(); P.set_real_precision(30)
CURVES = {'11a1': [0, -1, 1, -10, -20], '15a1': [1, 1, 1, -10, -10], '27a1': [0, 0, 1, 0, -7], '37a1': [0, 0, 1, -1, 0],
          '43a1': [0, 1, 1, 0, 0], '89a1': [1, 1, 1, -1, 0], '389a1': [0, 1, 1, -2, 0], 'c_A': [1, 0, 1, -3, 1],
          'c_B': [0, 0, 1, -3, 1], '5077a1': [0, 0, 1, -7, 6]}
OS = 32; AMAX = 3.0; P0 = 10_000; AMAXNUM = 1500

def spectrum(u, w, u0, u1, A, NT=12):
    dal = 1.0 / (OS * (u1 - u0)); h = 1.0 / (math.pi * A)
    M = 1 << int(math.ceil(math.log2(1.0 / (dal * h)))); h = 1.0 / (M * dal); J = int(A / dal)
    k = np.rint((u - u0) / h).astype(np.int64); d = (u - u0) - k * h
    out = np.zeros(J + 1, complex); al = np.arange(J + 1) * dal; dp = np.ones_like(d)
    for t in range(NT + 1):
        F = np.fft.fft(np.bincount(k, weights=w * dp, minlength=M))[:J + 1]
        out += ((-2j * math.pi * al) ** t / math.factorial(t)) * F; dp = dp * d
    return al, out * np.exp(-2j * math.pi * al * u0), dal

def truth(ai):
    E = P.ellinit(ai); N = int(P.ellglobalred(E)[0]); eps = int(P.ellrootno(E))
    L1 = float(P.real(P.lfun(P.lfunsympow(E, 2), 2))); bad = [int(x) for x in P.factor(N)[0]]; B = 1.0
    for p in P.primes(30000):
        p = int(p); Xp = 1.0 / p; ap = int(P.ellap(E, p))
        if p in bad:
            cp = p / (p - 1.0) if ap == 0 else 1 + 1.0 / p
            loc = float(P('my(E=ellinit(%s),f=lfuneuler(lfunsympow(E,2),%d)); subst(f,variable(f),1/%d^2)' % (ai, p, p)))
            B *= cp * loc
        else:
            l2 = ap * ap / p; B *= (1 - (l2 - 1) * Xp) / ((1 - Xp) * (1 - (l2 - 2) * Xp + Xp * Xp))
    return N, eps, B / L1

res = {}; t0 = time.time()
for name, ai in CURVES.items():
    if ONLY and name not in ONLY: continue
    v = P('my(E=ellinit(%s), Pr=primes([2,%d])); [Vecsmall(Pr), Vecsmall(apply(p->ellap(E,p), Pr))]' % (ai, X))
    p = np.array(list(v[0]), np.int64); a = np.array(list(v[1]), np.float64); lam = a / np.sqrt(p)
    small = {int(q): float(x) for q, x in zip(p[:2000], lam[:2000])}
    m = p >= P0; pf = p[m].astype(np.float64); u = np.sqrt(pf); u0, u1 = math.sqrt(P0), math.sqrt(X)
    w = lam[m] * np.log(pf) * S.hann(u, u0, u1); sig = math.sqrt(float(np.dot(w, w)) / 2)
    al, Sv, dal = spectrum(u, w, u0, u1, AMAX)
    Iu = np.linspace(u0, u1, 20001); unit = 1.5 * (2 * math.sqrt(2) / 3) * np.trapezoid(S.hann(Iu, u0, u1) * Iu ** 0.5, Iu)
    # good-prime factors g(l, ea, eb) relative to the r = 1 factor, from the data
    def lpow(l, k, bad):
        lp = small[l]
        if bad: return lp ** k
        t = [1.0, lp]
        for j in range(1, k): t.append(lp * t[j] - t[j - 1])
        return t[k]
    def factor_rel(l, ea, eb, bad, add):
        lp = small[l]
        if add:
            return 1.0 if (ea == 0 and eb == 0) else 0.0
        if bad:
            if eb == 0: return lpow(l, ea, True)              # A(k)/A(0) = lambda(l^k)
            if eb == 1: return (-lp / (l - 1.0)) / (1 + 1.0 / l)
            return 0.0
        rho = (lp * lp - 1.0) / l
        if eb == 0: Lv = lpow(l, ea, False) - lpow(l, ea + 2, False) / l
        elif eb == 1: Lv = -lp / l
        elif eb == 2: Lv = 1.0 / (l * (l - 1.0))
        else: Lv = 0.0
        return Lv / (1 - rho)
    def vv(x, l):
        k = 0
        while x % l == 0: x //= l; k += 1
        return k
    def fac(n):
        out = {}; q = 2
        while q * q <= n:
            while n % q == 0: out[q] = out.get(q, 0) + 1; n //= q
            q += 1
        if n > 1: out[n] = out.get(n, 0) + 1
        return out
    # all lines r = anum/bden (bden in 1..4, anum <= AMAXNUM), all-good weights and per-prime bad corrections
    lines = [(an, bd) for bd in (1, 2, 3, 4) for an in range(1, AMAXNUM + 1) if math.gcd(an, bd) == 1]
    rr = np.array([an / bd for an, bd in lines]); fg = np.zeros(len(lines)); facs = []
    for i, (an, bd) in enumerate(lines):
        ls = sorted(set(fac(an)) | set(fac(bd))); facs.append(ls); f = (an / bd) ** -0.25
        for l in ls: f *= factor_rel(l, vv(an, l), vv(bd, l), False, False)
        fg[i] = f
    corr = {}     # l -> (indices, bad-multiplicative weights, additive weights) for lines with l | an*bd
    for i, ls in enumerate(facs):
        for l in ls: corr.setdefault(l, []).append(i)
    corr_arr = {}
    for l, idx in corr.items():
        idx = np.array(idx); wb = np.zeros(len(idx)); wa = np.zeros(len(idx))
        for k, i in enumerate(idx):
            an, bd = lines[i]; f = (an / bd) ** -0.25; fa = f
            for l2 in facs[i]:
                gb = factor_rel(l2, vv(an, l2), vv(bd, l2), l2 == l, False)
                f *= gb; fa *= (factor_rel(l2, vv(an, l2), vv(bd, l2), True, True) if l2 == l else gb)
            wb[k] = f; wa[k] = fa
        corr_arr[l] = (idx, wb, wa)
    best = []
    ph = complex(math.cos(3 * math.pi / 4), math.sin(3 * math.pi / 4))
    for Np in range(2, NMAX + 1):
        badp = fac(Np); ws = fg.copy()
        for l in badp:
            if l in corr_arr:
                idx, wb, wa = corr_arr[l]
                ws[idx] = wa if abs(small.get(l, 0.0)) < 1e-12 else wb
                # lines divisible by two bad primes: recompute exactly (rare)
        multi = [l for l in badp if l in corr_arr]
        if len(multi) > 1:
            for l in multi:
                for i in corr_arr[l][0]:
                    an, bd = lines[i]
                    if sum(1 for l2 in facs[i] if l2 in badp) > 1:
                        f = (an / bd) ** -0.25
                        for l2 in facs[i]:
                            isb = l2 in badp; f *= factor_rel(l2, vv(an, l2), vv(bd, l2), isb, isb and abs(small.get(l2, 0.0)) < 1e-12)
                        ws[i] = f
        jj = np.rint(2 * np.sqrt(rr / Np) / dal).astype(np.int64)
        k = (jj < len(Sv)) & (ws != 0)
        if not k.any(): continue
        w2 = ws[k]; ss = Sv[jj[k]]
        proj = np.sum(w2 * ss) * ph; nrm2 = float(np.sum(w2 * w2))
        for eps in (1, -1):
            z = eps * proj.real / (sig * math.sqrt(nrm2))
            best.append((z, Np, eps, eps * proj.real / nrm2 / (unit * Np ** -0.25), int(k.sum())))
    best.sort(key=lambda x: -x[0])
    Nt, et, ct = truth(ai)
    ztrue = [b for b in best if b[1] == Nt and b[2] == et]
    rank = [i for i, b in enumerate(best) if b[1] == Nt and b[2] == et]
    res[name] = dict(N_true=Nt, eps_true=et, c_true=ct, N_hat=best[0][1], eps_hat=best[0][2], z_hat=best[0][0], c_hat=best[0][3],
                     n_lines=best[0][4], z_true=(ztrue[0][0] if ztrue else None), rank_true=(rank[0] if rank else None),
                     runners_up=[dict(N=b[1], eps=b[2], z=b[0]) for b in best[1:6]])
    print('%-7s true N=%5d eps=%+d c=%.3f | found N=%5d eps=%+d z=%.1f c_hat=%.3f (%d lines) | true at rank %s z=%s | next: %s  %.0fs'
          % (name, Nt, et, ct, best[0][1], best[0][2], best[0][0], best[0][3], best[0][4], res[name]['rank_true'],
             ('%.1f' % ztrue[0][0]) if ztrue else '-', ['%d%+d:%.1f' % (b[1], b[2], b[0]) for b in best[1:4]], time.time() - t0), flush=True)
    fo = os.path.join(S.HERE, 'out_spectrometer_%d.json' % X); old = json.load(open(fo)) if os.path.exists(fo) else {}
    old[name] = res[name]; json.dump(old, open(fo, 'w'), indent=1)
