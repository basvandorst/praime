"""Assemble the experiment's table from the saved runs (moment2.py / moment.py npz, mixed.py npz).
For each (N, k, T): diagonal D11 and its asymptotic, off-diagonal of M11 split into the arcs x = a c/q +- 1 (q | m) and the rest,
M2, L^2, X_Q, C2(Q), V(Q) = M11 + 2 L Re X_Q + L^2 C2(Q), V/T^2, V/D11; per-modulus arc weights against w(m)."""
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, glob, numpy as np
from smlib import Window, TWO_PI

def w_of_m(m):
    """w(m) = prod_{p^e || m} (1 + A_p/(p-1)^2), A_p = p-2 (e = 1), p-1 (e >= 2)."""
    r, x, p = 1.0, m, 2
    while p*p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0: x //= p; e += 1
            r *= 1 + ((p-2) if e == 1 else (p-1))/(p-1)**2
        p += 1
    if x > 1: r *= 1 + (x-2)/(x-1)**2
    return r

RUNS = {  # (N, k, T): (moment files, mixed file)
    (10007, 4, 6): (['mom_10007_4_6_primes.npz'], 'mixed_10007_4_6_primes.npz'),
    (10007, 4, 12): (['mom_10007_4_12_primes_1_100_3.5.npz'], 'mixed_10007_4_12_primes.npz'),
    (10007, 14, 12): (['mom_10007_14_12_primes_1_60_3.npz'], 'mixed_10007_14_12_primes.npz'),
    (10007, 14, 20): (['mom_10007_14_20_primes_1_150_3.npz'], 'mixed_10007_14_20_primes.npz'),
    (10007, 14, 30): (['mom_10007_14_30_primes_1_10_3.npz', 'mom_10007_14_30_primes_11_60_3.npz', 'mom_10007_14_30_primes_61_140_3.npz'], 'mixed_10007_14_30_primes.npz'),
    (5003, 14, 12): (['mom_5003_14_12_primes_1_80_3.npz'], 'mixed_5003_14_12_primes.npz'),
    (5003, 14, 20): (['mom_5003_14_20_primes_1_150_3.npz'], 'mixed_5003_14_20_primes.npz'),
    (5003, 14, 30): (['mom_5003_14_30_primes_1_12_3.npz', 'mom_5003_14_30_primes_13_150_3.npz'], 'mixed_5003_14_30_primes.npz'),
    (61, 4, 6): (['mom_61_4_6_primes_1_4000_3.5.npz'], 'mixed_61_4_6_primes.npz'),
}

def load(key):
    files, mixf = RUNS[key]
    rows = []; head = None
    for f in files:
        try: z = np.load(f, allow_pickle=True)
        except FileNotFoundError: return None
        if int(z['rows'][0, 0]) == 1: head = z
        rows.append(z['rows'])
    rows = np.concatenate(rows)
    try: mx = np.load(mixf)
    except FileNotFoundError: mx = None
    return head, rows, mx

def table(keys=None, verbose=True):
    out = {}
    for key in (keys or RUNS):
        r = load(key)
        if r is None: continue
        head, rows, mx = r
        N, k, T = key
        win = Window(N, k, T, h=0.05)
        diag_asym = float(np.sum(win.H**2*np.log(N*np.maximum(win.tt, 1e-9)**2/(4*math.pi**2)))*win.h/TWO_PI)
        b1, G1 = complex(head['b1']), complex(head['G1'])
        D11 = float(head['D11']); D2 = complex(head['D2'])
        off11 = rows[:, 3].sum(); off2 = complex(rows[:, 1].sum(), rows[:, 2].sum())
        M11 = D11 + off11 + abs(b1)**2 + 2*(np.conj(b1)*G1).real
        M2 = D2 + off2 + b1*b1 + 2*b1*G1
        has_arc = rows.shape[1] > 6
        arcs11 = rows[:, 8].sum() if has_arc else float('nan'); pm11 = rows[:, 11].sum() if has_arc else float('nan')
        L = float(head['L'])
        rec = dict(N=N, k=k, T=T, D11=D11, diag_asym=diag_asym, Dint=win.intH2/TWO_PI, off11=off11, arcs11=arcs11, pm11=pm11,
                   minor11=off11 - arcs11, M11=M11, M2=M2, L=L, Mmax=int(rows[:, 0].max()), rows=rows)
        if mx is not None:
            Qs = list(mx['res_Q']); i = Qs.index(10)
            XQ = complex(mx['res_XQ'][i]); C2 = float(mx['res_C2'][i])
            V = M11 + 2*L*XQ.real + L*L*C2
            rec.update(XQ=XQ, C2=C2, V=V)
        out[key] = rec
        if verbose:
            print('N=%5d k=%2d T=%2g (m<=%d): D11 = %.3f (asymptotic (1/2pi) int H^2 ln(Nt^2/4pi^2) = %.3f); off-diag M11 = %+.3f = arcs %+.3f '
                  '[x=+-1: %+.3f] + minor %+.3f; M11 = %.3f; M2 = %.3f%+.3fi; L^2 = %.3f; L^2 - (1/2pi)int H^2 = %.3f'
                  % (N, k, T, rec['Mmax'], D11, diag_asym, off11, arcs11, pm11, rec['minor11'], M11, M2.real, M2.imag, L*L, L*L - rec['Dint']))
            if 'V' in rec:
                print('        X_10 = %+.4f%+.4fi (-L C2 = %+.4f), C2(10) = %.4f, V(10) = %.3f, V/T^2 = %.4f, V/D11 = %.3f'
                      % (rec['XQ'].real, rec['XQ'].imag, -L*rec['C2'], rec['C2'], rec['V'], rec['V']/T**2, rec['V']/D11))
    return out


def lines_check(out):
    """X_q against the uniform Theorem D line sums -(1/2pi) sum_m mu(m) c_m(q^2)/(phi(m) q) H^(ln(m^2/q^2))."""
    from math import gcd
    def mob(x):
        r, m, p = 1, x, 2
        while p*p <= m:
            if m % p == 0:
                m //= p
                if m % p == 0: return 0
                r = -r
            p += 1
        return -r if m > 1 else r
    def ram(m, l): return sum(d*mob(m//d) for d in range(1, m+1) if m % d == 0 and l % d == 0)
    def phi(m): return sum(1 for a in range(1, m+1) if gcd(a, m) == 1)
    res = {}
    for key, rec in out.items():
        files, mixf = RUNS[key]
        try: mx = np.load(mixf)
        except FileNotFoundError: continue
        N, k, T = key; win = Window(N, k, T, h=0.05)
        Hh = lambda xi: np.sum(win.H*np.exp(-1j*xi*win.tt))*win.h
        worst = 0.0
        for q, xq in zip(mx['qs'], mx['Xq']):
            q = int(q)
            ln = -sum(mob(mm)*ram(mm, q*q)/(phi(mm)*q)*Hh(math.log(mm*mm/(q*q))) for mm in range(1, 41) if mob(mm) != 0)/TWO_PI
            worst = max(worst, abs(xq - ln))
        res[key] = worst
        print('N=%5d k=%2d T=%2g: max_q<=10 |X_q - Theorem D uniform line sum| = %.4f  (L = %.3f)' % (N, k, T, worst, rec['L']))
    return res

def per_modulus(out, mmax=4):
    for key, rec in out.items():
        rows = rec['rows']
        if rows.shape[1] <= 6: continue
        N, k, T = key; win = Window(N, k, T, h=0.05)
        s = []
        for m in range(1, mmax+1):
            r = rows[rows[:, 0] == m]
            if not len(r): continue
            r = r[0]; sh = int(round(TWO_PI*m/win.h))
            lag = np.sum(win.H[sh:]*win.H[:-sh])*win.h/math.pi if sh < len(win.H) else 0.0
            s.append('m=%d: %+.3f [arcs %+.3f, x=+-1 %+.3f, lag %+.3f]' % (m, r[3], r[8], r[11], lag))
        big = rows[(rows[:, 0] > mmax)]
        print('N=%5d k=%2d T=%2g M11 off-diagonal: ' % key + '; '.join(s) + '; m>%d: %+.3f (arcs %+.3f)' % (mmax, big[:, 3].sum(), big[:, 8].sum()))
        a2 = complex(rows[:, 6].sum(), rows[:, 7].sum()); p2 = complex(rows[:, 9].sum(), rows[:, 10].sum())
        wp = sum(w_of_m(int(r[0]))*complex(r[9], r[10]) for r in rows)
        print('        M2 off-diagonal %+.3f%+.3fi = arcs %+.3f%+.3fi + minor %+.3f%+.3fi ; x=+-1 parts %+.3f%+.3fi ; sum_m w(m)(x=+-1)_m = %+.3f%+.3fi'
              % (rows[:, 1].sum(), rows[:, 2].sum(), a2.real, a2.imag, rows[:, 1].sum() - a2.real, rows[:, 2].sum() - a2.imag, p2.real, p2.imag, wp.real, wp.imag))

if __name__ == '__main__':
    out = table()
    print()
    lines_check(out)
    print()
    per_modulus(out)
    print('\nPer-modulus arc weights (arcs / x=+-1 part) against w(m), where |x=+-1 part| > 0.05:')
    for key, rec in out.items():
        rows = rec['rows']
        if rows.shape[1] <= 6: continue
        for r in rows:
            m = int(r[0])
            for lab, a, p in (('M2', complex(r[6], r[7]), complex(r[9], r[10])), ('M11', r[8], r[11])):
                if abs(p) > 0.05:
                    print('  N=%d k=%d T=%g m=%d %s: arcs %.4f%+.4fi, x=+-1 %.4f%+.4fi, ratio %.3f%+.3fi, w(m) = %.4f'
                          % (key[0], key[1], key[2], m, lab, complex(a).real, complex(a).imag, complex(p).real, complex(p).imag,
                             (a/p).real, (a/p).imag, w_of_m(m)))
    M = 200000
    print('\nmean of w(m), m <= %d: %.5f ; zeta(2)/zeta(4) = 15/pi^2 = %.5f' % (M, np.mean([w_of_m(m) for m in range(1, M+1)]), 15/math.pi**2))
