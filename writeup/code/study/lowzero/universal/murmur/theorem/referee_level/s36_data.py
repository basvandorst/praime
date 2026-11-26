"""Referee checks S3 (normalisation) and S6 (own levels): window means of E_N(p) = trace(T_p W_N)-side / dim_new
from the trace formula (verified in s1_trace.py) with own Hurwitz class numbers (hcl.py), weight 2, p <= 4N, 16 windows.
Densities: own M^(N), prime-level M_2, Zubrilina's M_zub (her constants from the source)."""
import math, json, sys, time
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/referee_level')
from hcl import H_own, pari

PR = list(sp.primerange(2, 2 * 10 ** 6))
def u_(l): return 1.0 / ((l - 1) ** 2 * (l + 1))
C0 = math.exp(sum(math.log1p(-u_(l)) for l in PR))
def e_(l, r): return l * l / (l * l - 1.0) if r % l == 0 else 1.0 - u_(l)
def nup(r):
    v = C0
    for l in sp.factorint(r):
        v *= e_(l, 0) / (1 - u_(l))
    return v
def nuN(r, S):
    v = nup(r)
    for q in S:
        v /= e_(q, r)
    return v
def kap(S): return math.prod(q / (q - 1.0) for q in S)
def piN(S): return math.prod(1 + 1.0 / q for q in S)

# Zubrilina constants (source, Thm 1): products over p < 2e6 (tails O(1/X^2) or O(1/X))
la = np.array(PR, dtype=float)
alpha = 2 * math.pi * math.exp(np.sum(np.log((1 - la - 2 * la ** 2 + la ** 4) / (la ** 4 - 2 * la ** 2 + la))))
beta = 2 * math.pi * math.exp(np.sum(np.log((-1 + la ** 2 + la ** 3) / (la * (-1 + la + la ** 2)))))
gamma = 12 * math.exp(np.sum(np.log(la * (1 + la) / (-1 + la + la ** 2))))
def nuz(r): return math.prod(1 + l * l / (l ** 4 - 2 * l * l - l + 1.0) for l in sp.factorint(r)) if r > 1 else 1.0

def M_own(y, S):
    s = 2 * math.pi * piN(S) * math.sqrt(y) - 12 * kap(S) * y
    r = 1
    while r * r < 4 * y:
        s += (12 / math.pi) * kap(S) * nuN(r, S) * math.sqrt(4 * y - r * r)
        r += 1
    return s
def M_zub(y):
    s = beta * math.sqrt(y) - gamma * y
    r = 1
    while r * r < 4 * y:
        s += alpha * nuz(r) * math.sqrt(4 * y - r * r)
        r += 1
    return s
def wmean(f, a, b, n=4000):
    ys = a + (np.arange(n) + 0.5) * (b - a) / n
    return float(np.mean([f(y) for y in ys]))

def dim_new2(N):
    # own: g0(M) for squarefree M, Moebius-type inversion beta(q) = -2
    S = list(sp.factorint(N))
    def g0(M):
        P = list(sp.factorint(M)) if M > 1 else []
        mu = math.prod(q + 1 for q in P)
        n2 = math.prod(1 + (sp.legendre_symbol(-1 % q, q) if q > 2 else 0) for q in P)
        n3 = math.prod(1 + (sp.jacobi_symbol(-3, q) if q != 3 and q != 2 else (0 if q == 3 else -1)) for q in P)
        ninf = 2 ** len(P)
        return Fr(1) + Fr(mu, 12) - Fr(n2, 4) - Fr(n3, 3) - Fr(ninf, 2)
    tot = Fr(0)
    for d in sp.divisors(N):
        tot += (-2) ** len(sp.factorint(N // d)) * g0(d)
    return int(tot)

TYPES = {'2q': [2 * 10007, 2 * 10009, 2 * 10037], '3q': [3 * 6673, 3 * 6679, 3 * 6689],
         '6q': [6 * 3343, 6 * 3347, 6 * 3359], '5q': [5 * 4003, 5 * 4007, 5 * 4013],
         'S3': [2 * 5003, 6 * 1669]}
NW, HW = 16, 0.25

def level_data(N):
    S = list(sp.factorint(N))
    dim = dim_new2(N)
    dimp = int(pari(f'mfdim([{N},2],0)'))
    assert dim == dimp, (N, dim, dimp)
    fN = {int(q): 1 for q in S}
    num = np.zeros(NW); num0 = np.zeros(NW); den = np.zeros(NW)
    Ls = {r: [] for r in range(4)}
    for p in sp.primerange(2, 4 * N + 1):
        if N % p == 0:
            continue
        y = p / N
        w = min(int(math.ceil(y / HW)) - 1, NW - 1)
        tr = Fr(0); tr0 = None
        r = 0
        while r * r * N < 4 * p:
            m = 4 * p - r * r * N
            fm = pari.factor(m)
            fac = dict(fN)
            for i in range(len(fm[0])):
                l = int(fm[0][i]); fac[l] = fac.get(l, 0) + int(fm[1][i])
            D = N * m
            Hd = H_own(D, fac)
            Ls[r].append(math.pi * float(Hd) / math.sqrt(D))
            if r == 0:
                tr0 = Hd / 2; tr += Hd / 2
            else:
                tr += Hd
            r += 1
        tr -= p + 1
        lp = math.log(p)
        num[w] += lp * float(tr) / dim; num0[w] += lp * float(tr0) / dim; den[w] += lp
    return dict(N=N, S=[int(q) for q in S], dim=dim, data=list(num / den), data0=list(num0 / den),
                Lmean={r: (float(np.mean(v)), float(np.std(v) / math.sqrt(len(v))), len(v)) for r, v in Ls.items()})

t0 = time.time()
res = {}
for T, Ns in TYPES.items():
    for N in Ns:
        res[N] = level_data(N)
        print(f'level {N} {res[N]["S"]} dim {res[N]["dim"]} done {time.time()-t0:.0f}s', flush=True)
json.dump(dict(res=res, alpha=alpha, beta=beta, gamma=gamma, C0=C0),
          open('/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/referee_level/s36_data.json', 'w'))
print('alpha beta gamma C0', alpha, beta, gamma, C0)
