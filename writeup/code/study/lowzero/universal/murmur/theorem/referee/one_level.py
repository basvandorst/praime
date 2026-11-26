# Referee: one prime level N (not used in the note), weights 2 and 4, class-number form of sum_f eps_f sqrt(p) lambda_f(p)
# (T1, verified separately in t1_check.py), window means vs M_k (note) and Zubrilina's squarefree-level density.
import sys
import numpy as np
import cypari2
import sympy as sp

pari = cypari2.Pari()
pari.allocatemem(200 * 10**6)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 30011
YMAX, NW = 4.0, 16
h = YMAX / NW
assert pari.isprime(N)

# ---- Euler products (own code) ----
PR = np.array(list(sp.primerange(2, 2_000_000)), dtype=float)
def eprod(vals):
    return float(np.exp(np.sum(np.log(vals))))
q = PR
C0 = eprod(1 - 1 / ((q - 1) ** 2 * (q + 1)))         # prod (1 - u_l)
ALPHA = 2 * np.pi * eprod((1 - q - 2 * q**2 + q**4) / (q**4 - 2 * q**2 + q))
BETA = 2 * np.pi * eprod((-1 + q**2 + q**3) / (q * (-1 + q + q**2)))
GAMMA = 12 * eprod(q * (1 + q) / (-1 + q + q**2))

def nu_prime(r):
    v = C0
    for l in sp.factorint(r):
        v *= (l * l / (l * l - 1)) / (1 - 1 / ((l - 1) ** 2 * (l + 1)))
    return v

def nu_zub(r):
    v = 1.0
    for l in sp.factorint(r):
        v *= 1 + l**2 / (l**4 - 2 * l**2 - l + 1)
    return v

def U(n, x):
    # Chebyshev of the 2nd kind
    if n == 0:
        return np.ones_like(x)
    if n == 2:
        return 4 * x * x - 1
    raise ValueError

def M_note(y, k):
    y = np.asarray(y, float)
    s = np.zeros_like(y)
    for r in range(1, int(2 * np.sqrt(y.max())) + 2):
        arg = np.clip(4 * y - r * r, 0, None)
        s += nu_prime(r) * U(k - 2, r / (2 * np.sqrt(np.maximum(y, 1e-300)))) * np.sqrt(arg)
    out = (2 * np.pi * np.sqrt(y) + (-1) ** (k // 2 - 1) * (12 / np.pi) * s) / (k - 1)
    if k == 2:
        out -= 12 * y
    return out

def M_zub(y, k):
    y = np.asarray(y, float)
    s = np.zeros_like(y)
    for r in range(1, int(2 * np.sqrt(y.max())) + 2):
        arg = np.clip(4 * y - r * r, 0, None)
        s += nu_zub(r) * U(k - 2, r / (2 * np.sqrt(np.maximum(y, 1e-300)))) * np.sqrt(arg)
    out = (ALPHA * (-1) ** (k // 2 - 1) * s + BETA * np.sqrt(y)) / (k - 1)
    if k == 2:
        out -= GAMMA * y
    return out

def win_int(f, a, b, k, n=20001):
    x = np.linspace(a, b, n)
    v = f(x, k)
    return float(np.trapezoid(v, x) / (b - a))

# ---- data ----
dim = {2: int(pari(f"mfdim([{N},2],1)")), 4: int(pari(f"mfdim([{N},4],1)"))}
primes = [int(p) for p in pari(f"primes([2,{int(YMAX*N)}])") if int(p) != N]
E2, E4, ys, LR = [], [], [], {0: [], 1: [], 2: [], 3: []}
for p in primes:
    H0 = int(pari(f"qfbhclassno({4*p*N})"))
    s2 = H0 / 2 - (p + 1)
    s4 = H0 / 2
    LR[0].append((p, np.pi * H0 / np.sqrt(4 * p * N)))
    r = 1
    while r * r * N < 4 * p:
        D = 4 * p * N - r * r * N * N
        Hr = int(pari(f"qfbhclassno({D})"))
        x2 = r * r * N / (4 * p)
        s2 += Hr
        s4 += -(4 * x2 - 1) * Hr
        if r <= 3:
            LR[r].append((p, np.pi * Hr / np.sqrt(D)))
        r += 1
    E2.append(s2 / dim[2])
    E4.append(s4 / dim[4])
    ys.append(p / N)
ys = np.array(ys); w = np.log(np.array(primes, float)); E = {2: np.array(E2), 4: np.array(E4)}
print(f"N = {N}, dim S2 = {dim[2]}, dim S4 = {dim[4]}, primes p <= {YMAX:g}N: {len(primes)}, windows {NW} of h = {h}")

for k in (2, 4):
    chi_n, chi_z, rows = 0.0, 0.0, []
    for j in range(NW):
        a, b = j * h, (j + 1) * h
        m = (ys > a) & (ys <= b)
        ww, xx = w[m], E[k][m]
        mean = np.sum(ww * xx) / np.sum(ww)
        var = np.sum(ww**2 * (xx - mean) ** 2) / np.sum(ww) ** 2
        se = np.sqrt(var)
        pn, pz = win_int(M_note, a, b, k), win_int(M_zub, a, b, k)
        chi_n += ((mean - pn) / se) ** 2
        chi_z += ((mean - pz) / se) ** 2
        rows.append((a, b, int(m.sum()), mean, se, pn, pz))
    print(f"\nk = {k}: chi2 vs M_k (note) = {chi_n:.1f}, vs Zubrilina = {chi_z:.1f}  (on {NW} windows)")
    print("  window        n     mean       se      M_note    M_zub   (mean-M_note)/se")
    for a, b, n, mean, se, pn, pz in rows:
        print(f"  ({a:4.2f},{b:4.2f}] {n:5d} {mean:9.4f} {se:8.4f} {pn:9.4f} {pz:9.4f} {(mean-pn)/se:7.2f}")

print("\nmeans over primes of pi H(D_r)/sqrt(D_r):")
pred = {0: np.pi**2 / 6, 1: nu_prime(1), 2: nu_prime(2), 3: nu_prime(3)}
for r in range(4):
    v = np.array([x for _, x in LR[r]])
    z = (v.mean() - pred[r]) / (v.std() / np.sqrt(len(v)))
    print(f"  r = {r}: n = {len(v)}, mean {v.mean():.5f}, predicted {pred[r]:.5f}, z = {z:.2f}")
print(f"constants: C0 = nu'(1) = {C0:.9f}, alpha {ALPHA:.6f}, beta {BETA:.6f}, gamma {GAMMA:.6f}")
