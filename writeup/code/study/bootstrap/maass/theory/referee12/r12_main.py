# Referee code, stage 2 (own code).  Fourier inversion of ghat = Xi^2 B on stored GL nodes.
#   g(x) = (1/pi) int_0^inf ghat(t) cos(xt) dt
#   A^(k)(g) = (2/pi) int_0^inf ghat(t) [-log 2pi + Re psi(k/2+it)] dt
#   Delta(k) = A^(12) - A^(k) (t-route), and u-route 2 int_0^inf g(u) (e^{-ku/2}-e^{-6u})/(1-e^{-u}) du
import mpmath as mp, json, sys
from fractions import Fraction
mp.mp.dps = 60
HERE = '/home/bas/praime/writeup/code/study/bootstrap/maass/theory/'
sys.path.insert(0, HERE + 'referee12')
from r12_xi import TAU, gl_nodes

def load(NG):
    d = json.load(open(HERE + 'referee12/xi_nodes_%d.json' % NG))
    return [mp.mpf(v) for v in d['t']], [mp.mpf(v) for v in d['w']], [mp.mpf(v) for v in d['xi']]

def loadB(deg):
    d = json.load(open(HERE + 'cert12_d%d.json' % deg))
    return [mp.mpf(Fraction(s).numerator)/Fraction(s).denominator for s in d['b_exact']]

def primes(n):
    return [p for p in range(2, n+1) if all(p % q for q in range(2, int(p**0.5)+1))]

PR = primes(400)
NMAX = 400

class G:
    def __init__(self, NG, b):
        self.t, self.w, xi = load(NG)
        self.gh = [x*x*mp.fsum(bk*t**(2*k) for k, bk in enumerate(b)) for x, t in zip(xi, self.t)]
        self.cache = {}
    def g(self, x):
        x = mp.mpf(x)
        return mp.fsum(w*h*mp.cos(x*t) for w, h, t in zip(self.w, self.gh, self.t))/mp.pi
    def tint(self, f):
        return 2*mp.fsum(w*h*f(t) for w, h, t in zip(self.w, self.gh, self.t))/mp.pi
    def A(self, k):
        k = mp.mpf(k)
        return self.tint(lambda t: -mp.log(2*mp.pi) + mp.re(mp.digamma(k/2 + 1j*t)))
    def Delta(self, k):
        k = mp.mpf(k)
        return self.tint(lambda t: mp.re(mp.digamma(6 + 1j*t)) - mp.re(mp.digamma(k/2 + 1j*t)))
    def Aprime(self, s):  # (1/2pi) int_R ghat Re psi'(s/2+it) dt
        s = mp.mpf(s)
        return self.tint(lambda t: mp.re(mp.psi(1, s/2 + 1j*t)))/2
    def M(self, j):  # (1/2pi) int_R ghat t^{2j}
        return self.tint(lambda t: t**(2*j))/2

def cheb2(x, m):  # 2 T_m(x)
    c0, c1 = mp.mpf(2), 2*x
    if m == 0: return c0
    for _ in range(m-1):
        c0, c1 = c1, 2*x*c1 - c0
    return c1

def setup_P(G_):
    P = {}
    for p in PR:
        m, terms = 1, []
        while p**m <= NMAX:
            terms.append((m, mp.log(p)*mp.mpf(p)**(-mp.mpf(m)/2)*G_.g(m*mp.log(p))))
            m += 1
        P[p] = terms
    return P

def Pval(terms, x):
    return mp.fsum(c*cheb2(x, m) for m, c in terms)

def Pmin(terms):
    N = 2000
    xs = [-1 + mp.mpf(2)*i/N for i in range(N+1)]
    vs = [Pval(terms, x) for x in xs]
    cands = [(vs[0], xs[0]), (vs[-1], xs[-1])]
    for i in range(1, N):
        if vs[i] <= vs[i-1] and vs[i] <= vs[i+1]:
            a, b = xs[i-1], xs[i+1]
            for _ in range(160):  # golden section
                c1 = b - (b-a)/mp.phi; c2 = a + (b-a)/mp.phi
                if Pval(terms, c1) < Pval(terms, c2): b = c2
                else: a = c1
            xm = (a+b)/2; cands.append((Pval(terms, xm), xm))
    return min(cands)

def xp(p):
    return mp.mpf(TAU[p])/(2*mp.mpf(p)**(mp.mpf(11)/2))

def run(deg, NG):
    b = loadB(deg); Gb = G(NG, b)
    r = {'deg': deg, 'NG': NG}
    g0 = Gb.g(0); r['g0'] = mp.nstr(g0, 20)
    r['g(log n)/g0'] = {n: mp.nstr(Gb.g(mp.log(n))/g0, 6) for n in [2, 3, 4, 5, 7, 11, 19, 29, 97, 400]}
    P = setup_P(Gb)
    A12 = Gb.A(12)
    prim = 2*mp.fsum(Pval(P[p], xp(p)) for p in PR)
    r['A12/g0'] = mp.nstr(A12/g0, 20); r['Lemma1 residual/g0'] = mp.nstr((A12-prim)/g0, 5)
    # negative controls
    bad = 2*mp.fsum(Pval(P[p], xp(p) + (mp.mpf('0.001') if p == 2 else 0)) for p in PR)
    r['neg control x2+0.001 residual/g0'] = mp.nstr((A12-bad)/g0, 5)
    bad2 = 2*mp.fsum(Pval(P[p], (-xp(p) if p == 3 else xp(p))) for p in PR)
    r['neg control x3->-x3 residual/g0'] = mp.nstr((A12-bad2)/g0, 5)
    # Defect
    D, shares = 0, {}
    for p in PR:
        v, xm = Pmin(P[p]); dp = Pval(P[p], xp(p)) - v
        D += dp; shares[p] = (dp, xm)
    r['Defect/g0'] = mp.nstr(D/g0, 10)
    top = sorted(shares.items(), key=lambda kv: -kv[1][0])[:4]
    r['top shares (p, Defect_p/g0, argmin)'] = [(p, mp.nstr(v[0]/g0, 4), mp.nstr(v[1], 6)) for p, v in top]
    r['x_2, x_3, x_5'] = [mp.nstr(xp(p), 6) for p in (2, 3, 5)]
    r['min share negative?'] = mp.nstr(min(v[0] for v in shares.values())/g0, 3)
    # Delta and eps
    k = mp.findroot(lambda k: Gb.Delta(k) - 2*D, mp.mpf(12) - 2*D/Gb.Aprime(12))
    r['kappa* (Delta=2 Defect)'] = mp.nstr(k, 15); r['eps_exact'] = mp.nstr(12-k, 10)
    r["A'(12)/g0 exact"] = mp.nstr(Gb.Aprime(12)/g0, 10)
    # u-route at the root
    us, uw = gl_nodes(24)
    tot = 0
    for j in range(24):
        a0, a1 = mp.mpf(j)/4, mp.mpf(j+1)/4
        for x, w in zip(us, uw):
            u = a0 + (a1-a0)*(x+1)/2
            K = (mp.exp(-k*u/2) - mp.exp(-6*u))/(1 - mp.exp(-u))
            tot += w*(a1-a0)/2*Gb.g(u)*K
    r['u-route Delta(kappa*)/(2 Defect)'] = mp.nstr(2*tot/(2*D), 15)
    return r, Gb, D, g0

if __name__ == '__main__':
    deg = int(sys.argv[1]); NG = int(sys.argv[2])
    r, Gb, D, g0 = run(deg, NG)
    print(json.dumps(r, indent=1)); sys.stdout.flush()
