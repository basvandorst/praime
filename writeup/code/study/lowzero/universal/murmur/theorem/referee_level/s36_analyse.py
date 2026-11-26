import json, math, sys
import numpy as np
import sympy as sp
sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/referee_level')
D = json.load(open('/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/referee_level/s36_data.json'))
res = {int(k): v for k, v in D['res'].items()}
alpha, beta, gamma, C0 = D['alpha'], D['beta'], D['gamma'], D['C0']

def u_(l): return 1.0 / ((l - 1) ** 2 * (l + 1))
def e_(l, r): return l * l / (l * l - 1.0) if r % l == 0 else 1.0 - u_(l)
def nup(r):
    v = C0
    for l in sp.factorint(r):
        v *= e_(l, 0) / (1 - u_(l))
    return v
def nuN(r, S):
    v = nup(r)
    for q in S: v /= e_(q, r)
    return v
def kap(S): return math.prod(q / (q - 1.0) for q in S)
def piN(S): return math.prod(1 + 1.0 / q for q in S)
def nuz(r): return math.prod(1 + l * l / (l ** 4 - 2 * l * l - l + 1.0) for l in sp.factorint(r)) if r > 1 else 1.0
def M_own(y, S, k=2):
    s = 2 * math.pi * piN(S) * math.sqrt(y)
    r = 1
    sg = (-1) ** (k // 2 - 1)
    while r * r < 4 * y:
        x = r / (2 * math.sqrt(y)); U = 1.0 if k == 2 else 4 * x * x - 1
        s += sg * (12 / math.pi) * kap(S) * nuN(r, S) * U * math.sqrt(4 * y - r * r)
        r += 1
    return s / (k - 1) - (12 * kap(S) * y if k == 2 else 0)
def M_zub(y, k=2):
    s = beta * math.sqrt(y) / (k - 1) - (gamma * y if k == 2 else 0)
    r = 1; sg = (-1) ** (k // 2 - 1)
    while r * r < 4 * y:
        x = r / (2 * math.sqrt(y)); U = 1.0 if k == 2 else 4 * x * x - 1
        s += sg * alpha * nuz(r) * U * math.sqrt(4 * y - r * r) / (k - 1)
        r += 1
    return s
def wmean(f, a, b, n=3000):
    ys = a + (np.arange(n) + 0.5) * (b - a) / n
    return float(np.mean([f(y) for y in ys]))
NW, HW = 16, 0.25
W = [(w * HW, (w + 1) * HW) for w in range(NW)]
prime_pred = np.array([wmean(lambda y: M_own(y, []), a, b) for a, b in W])
zub_pred = np.array([wmean(M_zub, a, b) for a, b in W])

def own_pred(N):
    S = res[N]['S']; phi = math.prod(q - 1 for q in S)
    corr = (phi / 12) / res[N]['dim']
    return corr * np.array([wmean(lambda y: M_own(y, S), a, b) for a, b in W])

TYPES = {'2q': [20014, 20018, 20074], '3q': [20019, 20037, 20067], '6q': [20058, 20082, 20154], '5q': [20015, 20035, 20065]}
typeS = {'2q': [2], '3q': [3], '6q': [2, 3], '5q': [5]}
typepred = {T: np.array([wmean(lambda y: M_own(y, S), a, b) for a, b in W]) for T, S in typeS.items()}
# ---- between-level variances
v = {}; own = {}
for T, Ns in TYPES.items():
    dd = []
    for N in Ns:
        own[N] = own_pred(N)
        dd.append(np.array(res[N]['data']) - own[N])
    v[T] = np.var(np.array(dd), axis=0, ddof=1)
# pooled shape: v_{T,w} = c_T^2 s_w^2
s2 = np.ones(NW)
for it in range(20):
    c2 = {T: float(np.mean(v[T] / s2)) for T in TYPES}
    s2 = np.mean([v[T] / c2[T] for T in TYPES], axis=0)
    s2 /= np.mean(s2)
print('c_T =', {T: round(math.sqrt(c2[T]), 4) for T in TYPES})
print('s_w =', np.round(np.sqrt(s2), 3))
print('\nchi2 on 16 windows (sigma = c_T s_w from 3 levels of the type); rms in parentheses')
print('level      type | own        | prime-level | Zubrilina  | other types')
for T, Ns in TYPES.items():
    sig = np.sqrt(c2[T] * s2)
    for idx, N in enumerate(Ns):
        d = np.array(res[N]['data'])
        def chi(pred): return float(np.sum(((d - pred) / sig) ** 2)), float(np.sqrt(np.mean((d - pred) ** 2)))
        a, b, c = chi(own[N]), chi(prime_pred), chi(zub_pred)
        oth = ' '.join(f'{U}:{chi(typepred[U])[0]:.0f}' for U in TYPES if U != T)
        tag = 'TARGET' if idx == 0 else '      '
        print(f'{N} {tag} {T} | {a[0]:6.1f} ({a[1]:.3f}) | {b[0]:7.1f} ({b[1]:.3f}) | {c[0]:6.1f} ({c[1]:.3f}) | {oth}')
    dm = np.mean([res[N]['data'] for N in Ns], axis=0); om = np.mean([own[N] for N in Ns], axis=0)
    sm = sig / math.sqrt(3)
    print(f'   mean of 3 levels: own {np.sum(((dm-om)/sm)**2):.1f}, prime {np.sum(((dm-prime_pred)/sm)**2):.1f}, zub {np.sum(((dm-zub_pred)/sm)**2):.1f}')
# target levels: show data and predictions
for N in (20014, 20019, 20058, 20015):
    print(N, 'data', np.round(res[N]['data'], 3))
    print(N, 'own ', np.round(own[N], 3))
# ---- leave-one-type-out robustness: sigma from single-type variance only (no pooling)
print('\nunpooled sigma (per type, per window from 3 levels; 2 dof each): target chi2 own / prime / zub')
for T, Ns in TYPES.items():
    sig = np.sqrt(v[T]); N = Ns[0]; d = np.array(res[N]['data'])
    print(T, [round(float(np.sum(((d - p) / sig) ** 2)), 1) for p in (own[N], prime_pred, zub_pred)])

# ---- S3: normalisation: r = 0 part and class-number means
print('\nS3 normalisation: window means of (1/2)H(4pN)/dim vs 2 pi pi(N) sqrt(y) (and alternatives)')
for N in (10006, 10014, 20014, 20058, 20019, 20015):
    S = res[N]['S']; phi = math.prod(q - 1 for q in S); corr = (phi / 12) / res[N]['dim']
    d0 = np.array(res[N]['data0'])
    base = np.array([wmean(lambda y: 2 * math.pi * math.sqrt(y), a, b) for a, b in W]) * corr
    ratio = d0 / base
    print(f'  N={N} S={S}: fitted factor (data0 / 2 pi sqrt y) = {np.mean(ratio):.4f} +- {np.std(ratio)/4:.4f};'
          f' pi(N) = {piN(S):.4f}, kappa(N) = {kap(S):.4f}, 1')
    Lm = res[N]['Lmean']
    E0 = (math.pi ** 2 / 6) * math.prod(1 - q ** -2.0 for q in S)
    out = []
    for r in range(4):
        m, se, n = Lm[str(r)]
        pred = E0 if r == 0 else nuN(r, S)
        out.append(f'r={r}: {m:.4f}/{pred:.4f} (z {(m-pred)/se:+.2f}, n={n})')
    print('    mean pi H(D_r)/sqrt(D_r) / E^(N)(r):', '; '.join(out))
    own_ = own_pred(N); d = np.array(res[N]['data'])
    print(f'    whole: rms(data - own) = {np.sqrt(np.mean((d-own_)**2)):.3f}, rms(data - prime) = {np.sqrt(np.mean((d-prime_pred)**2)):.3f}, rms(data - zub) = {np.sqrt(np.mean((d-zub_pred)**2)):.3f}')
