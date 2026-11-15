"""Layout 2: why some moduli are lopsided (Rubinstein-Sarnak normalisation, Gaussian approximation).

Normalised race:  E(x) = (ln x / sqrt x) * h(x)/|R|  =  (r - 1)  +  sum_{chi real, nonprincipal mod P} sum_gamma x^(i gamma)/(1/2 + i gamma) + o(1)
(r = phi/|R| = number of real characters mod P).  Under GRH + LI the oscillating part has mean 0 and
variance  V(P) = sum_chi V(chi*),  V(chi*) = sum_gamma 1/(1/4 + gamma^2) over the zeros of the primitive L(s, chi*),
and (Davenport ch. 12, with Re(1/rho) = (1/2)/(1/4 + gamma^2) under GRH)
    V(chi_d) = ln(|d|/pi) + digamma((1 + a)/2) + 2 L'(1, chi_d)/L(1, chi_d),  a = 0 (d > 0), 1 (d < 0).
Gaussian approximation of the share of log-time with h < 0:  Phi(-(r - 1)/sqrt(V(P))).
The real characters mod P are the Kronecker symbols chi_d with d a fundamental discriminant, |d| dividing P.
Output: race_theory.json, race_theory.txt.
"""
import json
import math
import mpmath as mp
from scipy.stats import norm, spearmanr

mp.mp.dps = 40


def is_squarefree(n):
    n = abs(n)
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False
        p += 1
    return True


def is_fund(d):
    if d == 1 or d == 0:
        return False
    if d % 4 == 1:
        return is_squarefree(d)
    if d % 4 == 0:
        m = d // 4
        return m % 4 in (2, 3) and is_squarefree(m)
    return False


def kron(d, n):
    """Kronecker symbol (d/n), n >= 1."""
    res = 1
    while n % 2 == 0:
        n //= 2
        if d % 2 == 0:
            return 0
        res *= 1 if d % 8 in (1, 7) else -1
    # Jacobi symbol (d/n), n odd
    a = d % n
    t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return res * (t if n == 1 else 0)


fund = [d for d in range(-40, 41) if is_fund(d) and abs(d) <= 40]
V = {}
info = {}
for d in fund:
    q = abs(d)
    chi = [kron(d, n) if n > 0 else 0 for n in range(q)]
    # L(1) = -(1/q) sum chi(a) digamma(a/q)  (mpmath.dirichlet at s = 1 exactly raises the precision once per -1 value)
    L = -sum(chi[a] * mp.digamma(mp.mpf(a) / q) for a in range(1, q)) / q
    eps = mp.mpf(10) ** -12          # central difference (error ~ eps^2); mpmath's derivative at s = 1 is very slow
    Lp = (mp.dirichlet(1 + eps, chi) - mp.dirichlet(1 - eps, chi)) / (2 * eps)
    a = 0 if d > 0 else 1
    v = mp.log(mp.mpf(q) / mp.pi) + mp.digamma(mp.mpf(1 + a) / 2) + 2 * Lp / L
    V[d] = float(v)
    info[d] = {"L1": float(L), "Lp1": float(Lp), "V": float(v)}

out = []


def say(s=""):
    print(s)
    out.append(s)


say("Primitive real characters chi_d, |d| <= 40: V(chi_d) = sum over zeros 1/(1/4 + gamma^2)")
say("   d    L(1)     L'(1)    V")
for d in sorted(fund, key=abs):
    say(f"{d:4d} {info[d]['L1']:8.5f} {info[d]['Lp1']:8.5f} {info[d]['V']:8.5f}")
say()

main = json.load(open("race_main.json"))
rows = []
say(" P  r  characters (d)                 V(P)     z=(r-1)/sqrtV  Gauss P(h<0)  observed logdens<0 (10^7)  share primes<0")
for P in range(3, 41):
    ds = [d for d in fund if P % abs(d) == 0]
    r = main[str(P)]["r"]
    assert len(ds) == r - 1, (P, ds, r)
    VP = sum(V[d] for d in ds)
    z = (r - 1) / math.sqrt(VP)
    g = norm.cdf(-z)
    o = main[str(P)]["10000000"]
    rows.append((P, r, ds, VP, z, g, o["logdens_neg"], o["share_primes_neg"]))
    say(f"{P:2d} {r:2d}  {str(ds):30s} {VP:7.4f}  {z:6.3f}        {g:.2e}      {o['logdens_neg']:.5f}              {o['share_primes_neg']:.5f}")
rho = spearmanr([x[4] for x in rows], [-x[6] for x in rows])
say()
say(f"Spearman rank correlation between z and -(observed log density of h < 0): {rho.correlation:.3f} (p = {rho.pvalue:.1e})")
never = [x[0] for x in rows if main[str(x[0])]["10000000"]["first_neg"] is None]
some = [x[0] for x in rows if main[str(x[0])]["10000000"]["first_neg"] is not None]
say(f"min z among never-negative P: {min(x[4] for x in rows if x[0] in never):.3f}; "
    f"max z among P that go negative: {max(x[4] for x in rows if x[0] in some):.3f}")
say("P sorted by z: " + ", ".join(f"{x[0]}:{x[4]:.2f}" for x in sorted(rows, key=lambda x: -x[4])))
json.dump({"V": {str(d): info[d] for d in fund},
           "rows": [{"P": x[0], "r": x[1], "ds": x[2], "V": x[3], "z": x[4], "gauss_neg": x[5]} for x in rows]},
          open("race_theory.json", "w"), indent=1)
open("race_theory.txt", "w").write("\n".join(out) + "\n")
