# Own recomputation of the one-point sums from the zero files (no repo code imported).
# theta_L(t) = (t/2) ln(Q/pi^d) + sum_j Im lnGamma((1/2 + mu_j + i t)/2)
import os, sys, json, numpy as np
from scipy.special import loggamma
D = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/'
OUT = os.path.dirname(os.path.abspath(__file__))+'/'

# my own data: (Q, mu, eps) derived independently; eps of complex characters from Gauss sums below
def gauss_eps(q, gen_val):
    # character mod prime q defined by chi(g) = gen_val for primitive root g
    import sympy
    g = sympy.primitive_root(q)
    chi = {}
    x = 1
    for k in range(q - 1):
        chi[x] = gen_val ** k
        x = x * g % q
    tau = sum(chi[a] * np.exp(2j * np.pi * a / q) for a in range(1, q))
    a = 0 if abs(chi[q - 1] - 1) < 1e-9 else 1
    return chi, tau / ((1j) ** a * np.sqrt(q)), a

def theta(t, Q, mu):
    t = np.asarray(t, dtype=float)
    d = len(mu)
    s = 0.5 * t * np.log(Q / np.pi ** d)
    for m in mu:
        s = s + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return s

def theta_p(t, Q, mu, h=1e-5):
    return (theta(t + h, Q, mu) - theta(t - h, Q, mu)) / (2 * h)

def own_data():
    L = {}
    for k in (12, 16, 18, 20, 22, 26):
        key = 'D12' if k == 12 else 'F%d' % k
        L[key] = dict(Q=1.0, mu=[(k - 1) / 2, (k + 1) / 2], eps=complex((1j) ** k), fam='GL2 level 1')
    for key, N in (('R23', 23), ('R31', 31)):
        L[key] = dict(Q=float(N), mu=[0.0, 1.0], eps=1 + 0j, fam='GL2 Artin S3')
    for key, N in (('S11a', 11), ('S14a', 14), ('S15a', 15), ('S17a', 17), ('S19a', 19), ('S37a', 37)):
        L[key] = dict(Q=float(N * N), mu=[1.0, 1.0, 2.0], eps=1 + 0j, fam='GL3 Sym2 E')
    for key, N in (('G249', 249), ('G277', 277)):
        L[key] = dict(Q=float(N), mu=[0.5, 0.5, 1.5, 1.5], eps=1 + 0j, fam='GL4 genus 2')
    # characters: X5 = Conrey 5.2 (chi(2) = i), X5d = 5.3 (chi(2) = -i); X7 = Conrey 7.2, X7d = 7.4
    # Conrey 7.2: chi(3) = e(2/6)?? -- determine by matching: Conrey label q.n: chi(g^k) for g primitive root ...
    return L

if __name__ == '__main__':
    pred = json.load(open(D + 'predictions.json'))['Lfunctions']
    L = own_data()
    # complex characters: compute eps for both choices, match to frozen value only to identify which is which
    for key in ('X5', 'X5d', 'X7', 'X7d'):
        q = 5 if key.startswith('X5') else 7
        P = pred[key]
        cands = []
        import cmath
        order = 4 if q == 5 else 3
        for j in range(1, order):
            if np.gcd(j, order) != 1:
                continue
            chi, e, a = gauss_eps(q, cmath.exp(2j * cmath.pi * j / (q - 1) * ((q - 1) // order)))
            cands.append((e, a, j))
        fe = complex(*P['eps'])
        best = min(cands, key=lambda c: abs(c[0] - fe))
        L[key] = dict(Q=float(q), mu=[float(best[1])], eps=best[0], fam='GL1 complex character')
        print(key, 'own eps candidates', [(np.round(c[0], 6), c[1]) for c in cands], 'frozen', P['eps'], 'mu frozen', P['mu'])
    res = {}
    for key, v in L.items():
        P = pred[key]
        assert abs(P['Q'] - v['Q']) < 1e-9 and np.allclose(P['mu'], v['mu']), (key, P['Q'], P['mu'])
        assert abs(complex(*P['eps']) - v['eps']) < 1e-6, (key, P['eps'], v['eps'])
        g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']
        T2 = float(P['T2'])
        g = np.sort(g[(g > 0) & (g <= T2)])
        th = theta(g, v['Q'], v['mu'])
        G = np.sum(np.exp(2j * th))
        c = P['c']
        pr = -v['eps'] * c * T2 / (2 * np.pi)
        n = len(g)
        # counting: n - theta(T2)/pi; offsets
        off = np.arange(1, n + 1) - theta(g, v['Q'], v['mu']) / np.pi
        res[key] = dict(n=n, G=[G.real, G.imag], pred=[pr.real, pr.imag], slope=float((G / pr).real),
                        end=float(n - theta(T2, v['Q'], v['mu']) / np.pi),
                        off_mean=float(off.mean()), drift=float(off[-n // 7:].mean() - off[:n // 7].mean()))
        print('%-5s n=%4d G=%9.2f%+9.2fi pred=%9.2f%+9.2fi slope=%.4f  n-theta(T2)/pi=%+.2f  mean(n-theta(g_n)/pi)=%.3f drift=%+.3f'
              % (key, n, G.real, G.imag, pr.real, pr.imag, (G / pr).real, res[key]['end'], res[key]['off_mean'], res[key]['drift']))
    json.dump(res, open(OUT + 'r1_sums.json', 'w'), indent=1)
