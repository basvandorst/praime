"""theory.py -- the linearised explicit formula for prime patterns.

Pattern: linear forms L_i(n) = a_i n + b_i (i = 1..k).  E(x) = sum_{n<=x} prod_i Lambda(L_i(n)) - S x.
Conjecture (linear terms in the zeros):
    E(x) ~ S * sum_i (1/a_i) * [ (psi(a_i x) - a_i x) + sum_{chi* != 1 primitive} c_i(chi*) psi(a_i x, chi*) ]
with c_i(chi*) = mean of conj(chi*(L_i(n))) over the n mod q* (q* = conductor) with gcd(prod_j L_j(n), q*) = 1.
Each zero rho of L(s, chi*) then enters E(x) as  - A(rho) x^rho / rho  with  A(rho) = S * sum_i c_i(chi*) a_i^(rho - 1).
"""
import numpy as np
from math import gcd
from fractions import Fraction
import sympy

C2 = 0.66016181584686957392781211001455577843262336028473


def singular_series(forms, P=2_000_000):
    """Hardy-Littlewood constant S = prod_p (1 - nu_p/p) / (1 - 1/p)^k for the forms (a_i, b_i)."""
    k = len(forms)
    primes = list(sympy.primerange(2, P))
    logS = 0.0
    for p in primes:
        if p < 200 or any(a % p == 0 for a, _ in forms):
            nu = sum(1 for n in range(p) if any((a * n + b) % p == 0 for a, b in forms))
        else:
            roots = set()
            for a, b in forms:
                roots.add((-b * pow(a, -1, p)) % p)
            nu = len(roots)
        logS += np.log1p(-nu / p) - k * np.log1p(-1.0 / p)
    # tail for p >= P: nu = k, log factor ~ -binom(k,2)/p^2 ; sum_{p>P} 1/p^2 ~ 1/(P ln P)
    logS += -(k * (k - 1) / 2) / (P * np.log(P))
    return float(np.exp(logS))


# ---------------- Dirichlet characters ----------------
def _prime_power_chars(p, e):
    q = p ** e
    units = [n for n in range(q) if n % p]
    out = []
    if p == 2:
        if e == 1:
            return [np.array([0, 1], dtype=complex)]
        # n = (-1)^a 5^b mod 2^e
        m = 1 << max(e - 2, 0)
        ab = {}
        for a in (0, 1):
            for b in range(m):
                ab[(pow(-1, a) * pow(5, b, q)) % q] = (a, b)
        for s in (0, 1):
            for t in range(m):
                v = np.zeros(q, dtype=complex)
                for n in units:
                    a, b = ab[n]
                    v[n] = (-1) ** (a * s) * np.exp(2j * np.pi * t * b / m)
                out.append(v)
        return out
    g = sympy.primitive_root(q)
    phi = q - q // p
    dlog = {}
    x = 1
    for i in range(phi):
        dlog[x] = i
        x = x * g % q
    for k in range(phi):
        v = np.zeros(q, dtype=complex)
        for n in units:
            v[n] = np.exp(2j * np.pi * k * dlog[n] / phi)
        out.append(v)
    return out


def characters(q):
    """All Dirichlet characters mod q as value arrays of length q."""
    fac = sympy.factorint(q)
    comps = [(p ** e, _prime_power_chars(p, e)) for p, e in fac.items()]
    chars = [np.ones(q, dtype=complex)]
    mods = [1]
    for pe, cl in comps:
        new = []
        for v0 in chars:
            for v1 in cl:
                v = np.array([v0[n % len(v0)] * v1[n % pe] for n in range(q)]) if len(v0) > 1 else np.array([v1[n % pe] for n in range(q)])
                new.append(v)
        chars = new
    for v in chars:
        for n in range(q):
            if gcd(n, q) != 1:
                v[n] = 0
    return chars


def conductor(v):
    q = len(v)
    for d in sorted(sympy.divisors(q)):
        if all(abs(v[n] - 1) < 1e-9 for n in range(q) if gcd(n, q) == 1 and n % d == 1 % d):
            return d
    return q


def primitive_characters(qmax=30):
    """List of dicts: q (conductor), v (values mod q), parity, order, real, label. One per character (both of a conjugate pair kept)."""
    out = []
    for q in range(3, qmax + 1):
        for v in characters(q):
            if conductor(v) != q:
                continue
            vals = v[[n for n in range(q) if gcd(n, q) == 1]]
            order = 1
            while not np.allclose(vals ** order, 1):
                order += 1
            par = int(round(v[q - 1].real))
            real = bool(np.allclose(v.imag, 0))
            out.append(dict(q=q, v=np.round(v, 12), parity=par, order=order, real=real))
    # labels: q.index ; mark conjugate pairs
    byq = {}
    for c in out:
        byq.setdefault(c["q"], []).append(c)
    for q, cl in byq.items():
        # sort: real first, then by order, then by argument of value at the smallest generator-like unit
        g = next(n for n in range(2, q) if gcd(n, q) == 1 and all(not np.isclose(c["v"][n], 1) or c["order"] == 1 for c in cl) ) if q > 2 else 1
        cl.sort(key=lambda c: (not c["real"], c["order"], np.angle(c["v"][2 if gcd(2, q) == 1 else 7 % q]) % (2 * np.pi)))
        for i, c in enumerate(cl):
            c["label"] = "%d.%d" % (q, i + 1)
    for c in out:
        cc = np.conj(c["v"])
        for d in byq[c["q"]]:
            if np.allclose(d["v"], cc):
                c["conj"] = d["label"]
    return out


def local_coeffs(forms, chi):
    """c_i(chi) = mean over admissible n mod q of conj(chi(a_i n + b_i))."""
    q = len(chi["v"])
    adm = [n for n in range(q) if gcd(int(np.prod([(a * n + b) % q for a, b in forms])) % q if q > 1 else 1, q) == 1]
    if not adm:
        return [0j] * len(forms)
    return [complex(np.mean([np.conj(chi["v"][(a * n + b) % q]) for n in adm])) for a, b in forms]


PATTERNS = {
    "h2": [(1, 0), (1, 2)],
    "h4": [(1, 0), (1, 4)],
    "h6": [(1, 0), (1, 6)],
    "h8": [(1, 0), (1, 8)],
    "h10": [(1, 0), (1, 10)],
    "h12": [(1, 0), (1, 12)],
    "h30": [(1, 0), (1, 30)],
    "SG": [(1, 0), (2, 1)],
    "T026": [(1, 0), (1, 2), (1, 6)],
    "T046": [(1, 0), (1, 4), (1, 6)],
}


def singular_const(name):
    if name.startswith("h"):
        h = int(name[1:])
        s = 2 * C2
        for p in sympy.primefactors(h):
            if p > 2:
                s *= (p - 1) / (p - 2)
        return s
    if name == "SG":
        return 2 * C2
    return singular_series(PATTERNS[name])


def scales(name):
    return sorted(set(a for a, b in PATTERNS[name]))


def prediction_table(name, chars):
    """For each primitive character and each scale a (= the a_i of the forms): the coefficient
    beta_a(chi) = (S/a) * sum_{i: a_i = a} c_i(chi) on psi(a x, chi). (Shifts b_i only move psi(a x + b_i) by O(b_i ln x).)"""
    forms = PATTERNS[name]
    S = singular_const(name)
    rows = []
    for c in chars:
        ci = local_coeffs(forms, c)
        beta = {a: S / a * sum(x for x, (a2, b) in zip(ci, forms) if a2 == a) for a in scales(name)}
        beta = {a: (0j if abs(v) < 1e-12 else v) for a, v in beta.items()}
        rows.append(dict(label=c["label"], q=c["q"], parity=c["parity"], order=c["order"], real=c["real"],
                         conj=c.get("conj"), c=ci, beta=beta, present=any(abs(v) > 1e-12 for v in beta.values())))
    return S, rows


if __name__ == "__main__":
    chars = primitive_characters(30)
    print("primitive characters of conductor <= 30:", len(chars))
    for name in PATTERNS:
        S, rows = prediction_table(name, chars)
        forms = PATTERNS[name]
        zeta = {a: S / a * sum(1 for a2, b in forms if a2 == a) for a in scales(name)}
        print("\n== pattern %s forms %s  S = %.10f  zeta: coefficient on psi(a x) - a x: %s" % (name, forms, S, {a: round(v, 6) for a, v in zeta.items()}))
        for r in rows:
            if r["present"]:
                bs = "  ".join("a=%d: %+.4f%+.4fi (rel %.4f)" % (a, v.real, v.imag, abs(v) / zeta[a]) for a, v in r["beta"].items())
                print("  %-6s q=%2d %s order %2d %s  beta: %s" % (r["label"], r["q"], "even" if r["parity"] == 1 else "odd ", r["order"], "real" if r["real"] else "cplx", bs))
        absent = [r["label"] for r in rows if not r["present"]]
        print("  absent (%d): %s" % (len(absent), " ".join(absent)))
