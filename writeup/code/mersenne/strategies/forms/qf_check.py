"""Independent check of qf_count.py (questions 1 and 2), method B.

- primes from sympy.primerange (not the numpy sieve),
- representability of each prime p as a^2 + P b^2 by Cornacchia's algorithm
  (square root of -P mod p via sympy.sqrt_mod, then the Euclidean descent),
  with a direct verification x^2 + P y^2 == p of every solution found,
- h(-4P) and the class group structure from PARI (cypari2: qfbclassno and
  quadclassunit), not from counting reduced forms,
- idoneal test by genus theory: P idoneal iff h(-4P) equals the number of
  genera 2^(mu-1) (Cox, Primes of the form x^2+ny^2, Prop. 3.11 / Thm 3.15),
  and also iff the class group has exponent <= 2 (PARI cyc),
- the residue-class test redone from the Cornacchia results.
Compares everything with qf_count.json and prints the number of mismatches.
"""
import json
from math import gcd, isqrt

import cypari2
from sympy import primerange, primefactors
from sympy.ntheory import sqrt_mod

pari = cypari2.Pari()
N = 500_000


def cornacchia_prime(P, p):
    """Return (x, y) with x^2 + P y^2 = p, or None.  p prime."""
    if p == 2:
        for y in range(2):
            x2 = 2 - P * y * y
            if x2 >= 0 and isqrt(x2) ** 2 == x2:
                return (isqrt(x2), y)
        return None
    if P % p == 0:
        return (0, 1) if P == p else None
    if P > p:
        return None
    if pow((-P) % p, (p - 1) // 2, p) != 1:
        return None
    r0 = sqrt_mod((-P) % p, p)
    for r in (r0, p - r0):
        a, b = p, r
        lim = isqrt(p)
        while b > lim:
            a, b = b, a % b
        rest = p - b * b
        if rest % P == 0:
            c = rest // P
            y = isqrt(c)
            if y * y == c:
                assert b * b + P * y * y == p
                return (b, y)
    return None


def genera(P):
    """Number of genera of forms of discriminant -4P (Cox Prop. 3.11)."""
    r = len([q for q in primefactors(P) if q != 2])
    if P % 4 == 3:
        mu = r
    elif P % 4 in (1, 2):
        mu = r + 1
    elif P % 8 == 4:
        mu = r + 1
    else:  # P % 8 == 0
        mu = r + 2
    return 2 ** (mu - 1)


def main():
    ref = json.load(open("qf_count.json"))
    primes = list(primerange(2, N + 1))
    assert len(primes) == ref["pi_N"] == 41538
    bad = 0
    out = []
    for row in ref["rows"]:
        P = row["P"]
        rep = {}
        for p in primes:
            rep[p] = cornacchia_prime(P, p) is not None
        drawn_all = sum(rep.values())
        excl = [p for p in primes if (4 * P) % p != 0]
        drawn_excl = sum(rep[p] for p in excl)
        h_pari = int(pari.qfbclassno(-4 * P))
        cyc = [int(c) for c in pari.quadclassunit(-4 * P)[1]]
        exp2 = all(c <= 2 for c in cyc)
        g = genera(P)
        idoneal_genus = (h_pari == g)
        m = 4 * P
        state = {}
        for p in excl:
            s = state.setdefault(p % m, set())
            s.add(rep[p])
        decided = all(len(s) == 1 for s in state.values()) and \
            len(state) == sum(1 for x in range(m) if gcd(x, m) == 1)
        checks = {
            "drawn_all": drawn_all == row["drawn_all"],
            "drawn_excl": drawn_excl == row["drawn_excl"],
            "h": h_pari == row["h"],
            "decided": decided == row["decided"],
            "idoneal(genus)": idoneal_genus == row["idoneal_list"],
            "idoneal(exp2)": exp2 == row["idoneal_list"],
            "decided==idoneal": decided == idoneal_genus,
        }
        nb = sum(not v for v in checks.values())
        bad += nb
        out.append({"P": P, "h_pari": h_pari, "cyc": cyc, "genera": g,
                    "drawn_all": drawn_all, "decided": decided,
                    "idoneal_genus": idoneal_genus})
        print(f"P={P:3d} h={h_pari:2d} cyc={cyc} genera={g} drawn={drawn_all} "
              f"decided={decided} idoneal={idoneal_genus} "
              f"{'OK' if nb == 0 else 'MISMATCH ' + str(checks)}", flush=True)
    print("total mismatches:", bad)
    json.dump(out, open("qf_check.json", "w"))


if __name__ == "__main__":
    main()
