"""Recheck of bias.py from qf_check.json only (Cornacchia counts, PARI h,
genus-theory genera count), primes up to sqrt N from sympy.
drawn_excl = drawn_all - (2 drawn ? 1 : 0) - (odd prime P drawn ? 1 : 0):
2 = x^2 + P y^2 only for P = 1, 2; an odd prime q | P is drawn only when q = P.
"""
import json
from math import isqrt

from sympy import isprime, primepi

N = 500_000


def main():
    ref = json.load(open("bias.json"))["rows"]
    rows = json.load(open("qf_check.json"))
    pis = int(primepi(isqrt(N)))
    n_all = int(primepi(N))
    bad = 0
    tot_o = tot_p = 0.0
    pos = 0
    for r, f in zip(rows, ref):
        P, h, g = r["P"], r["h_pari"], r["genera"]
        special = [q for q in range(2, P + 1) if (4 * P) % q == 0 and isprime(q)]
        if 2 not in special:
            special.append(2)
        drawn_special = (1 if P in (1, 2) else 0) + (1 if isprime(P) and P > 2 else 0)
        n_excl = n_all - len(set(special))
        drawn_excl = r["drawn_all"] - drawn_special
        obs = n_excl / (2 * h) - drawn_excl
        pred = 0.5 * (g + h - 1) / (2 * h) * pis
        bad += abs(obs - f["observed_deficit"]) > 1e-9 or abs(pred - f["predicted_mean_deficit"]) > 1e-9
        tot_o += obs
        tot_p += pred
        pos += obs > 0
    print({"mismatches": bad, "P with deficit > 0": pos, "sum observed": tot_o,
           "sum predicted": tot_p})


if __name__ == "__main__":
    main()
