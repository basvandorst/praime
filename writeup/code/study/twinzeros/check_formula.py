"""check_formula.py -- three independent computations of the coefficient of x^rho/rho for pairs (n, n+h):
 (1) local densities (theory.local_coeffs, enumeration of admissible classes),
 (2) the closed form  beta = 2 S(h) mu(q*) conj(chi(h)) / prod_{p | q*} (p - 2)  for even chi of odd squarefree conductor
     coprime to h, 0 otherwise,
 (3) the circle method: major arcs a/q, Gauss sums of the induced characters, summed over q = q* m (q <= QMAX).
"""
import numpy as np, sympy
from math import gcd
import theory

chars = theory.primitive_characters(30)
byl = {c["label"]: c for c in chars}


def closed(c, h, S):
    q = c["q"]
    if c["parity"] != 1 or q % 2 == 0 or not sympy.factorint(q) or max(sympy.factorint(q).values()) > 1 or gcd(q, h) > 1:
        return 0j
    pr = np.prod([p - 2 for p in sympy.primefactors(q)])
    return 2 * S * sympy.mobius(q) * np.conj(c["v"][h % q]) / pr


# (1) vs (2) for all 167 characters and seven gaps
worst = 0.0
for h in (2, 4, 6, 8, 10, 12, 30):
    S, rows = theory.prediction_table("h%d" % h, chars)
    for r in rows:
        d = abs(r["beta"][1] - closed(byl[r["label"]], h, S))
        worst = max(worst, d)
print("local densities vs closed form, 167 characters x 7 gaps: max |difference| = %.2e" % worst)

# (3) circle method
QMAX = 900
HS3 = (2, 4, 6)
def induced(c, q):
    n = np.arange(q)
    v = c["v"][n % c["q"]].copy()
    v[np.gcd(n, q) != 1] = 0
    return v

def coef_circle(c):
    """coefficient of N^rho/rho (rho a zero of L(s, chi)) in sum Lambda(n)Lambda(n+h), h in HS3, from the major arcs q <= QMAX."""
    qs = c["q"]
    cc = byl[c["conj"]] if not c["real"] else c
    tot = np.zeros(len(HS3), dtype=complex)
    for q in range(qs, QMAX + 1, qs):
        mu = int(sympy.mobius(q))
        if mu == 0:
            continue
        n = np.arange(q)
        units = n[np.gcd(n, q) == 1]
        phi = len(units)
        E = np.exp(2j * np.pi * np.outer(units, n) / q)        # e(a b / q)
        c_chi = (E @ np.conj(induced(c, q))) / phi             # c_chi(a) = (1/phi) sum_b conj(chi(b)) e(ab/q)
        c_chib = (E @ np.conj(induced(cc, q))) / phi
        main = mu / phi
        for i, h in enumerate(HS3):
            eh = np.exp(-2j * np.pi * h * units / q)
            # fluct_chi x conj(main) gives N^rho/rho; main x conj(fluct of the conjugate character) gives N^rho/rho too
            tot[i] += -np.sum(c_chi * main * eh) - np.sum(main * np.conj(c_chib) * eh)
    return tot

print("circle method (q <= %d) vs local densities: coefficient of x^rho/rho, which should equal -beta" % QMAX)
tabs = {h: {r["label"]: r for r in theory.prediction_table("h%d" % h, chars)[1]} for h in HS3}
for lab in ("3.1", "4.1", "5.1", "5.2", "7.1", "7.2", "7.4", "15.1", "15.2", "21.1", "21.3", "8.1", "12.1", "13.2"):
    c = byl[lab]
    cm = coef_circle(c)
    for i, h in enumerate(HS3):
        beta = tabs[h][lab]["beta"][1]
        print("  %-5s %s h=%d  circle %+.4f%+.4fi   -beta(local) %+.4f%+.4fi" % (lab, "even" if c["parity"] == 1 else "odd ", h, cm[i].real, cm[i].imag, -beta.real, -beta.imag))
