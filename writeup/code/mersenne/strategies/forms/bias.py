"""Why every P has fewer drawn primes than 1/(2h) predicts (Chebyshev-type bias).

Heuristic (prime-square term of the explicit formula, as in Rubinstein-Sarnak;
for Chebotarev sets see Ng 2000 / Bellaiche 2016): the primes of the principal
form are those whose Frobenius in G = Gal(L/Q) (L the ring class field, G
generalized dihedral of order 2h) is the identity.  Counting with weights 1/k
over prime powers is unbiased; removing the prime squares p^2 <= x, whose
Frobenius squared is 1 for a share I/|G| of p, gives
    pi_1(x) - pi(x)/(2h)  ~  -(1/2) * (I - 1)/(2h) * pi(sqrt x),
I = #{g in G: g^2 = 1} = #Cl[2] + h  (Cl[2] = number of genera, h reflections).
Compares this mean prediction with the observed deficit (n_excl/(2h) - drawn_excl).
Inputs: qf_count.json (counts, h), PARI quadclassunit for #Cl[2].
This is a heuristic mean; the actual value at one x fluctuates (zeros).
"""
import json
from math import isqrt

import cypari2
import numpy as np

pari = cypari2.Pari()
N = 500_000


def main():
    rows = json.load(open("qf_count.json"))["rows"]
    s = np.ones(isqrt(N) + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(isqrt(N)) + 1):
        if s[i]:
            s[i * i::i] = False
    pi_sqrt = int(s.sum())
    out = []
    obs_all, pred_all = [], []
    for r in rows:
        P, h = r["P"], r["h"]
        cyc = [int(c) for c in pari.quadclassunit(-4 * P)[1]]
        cl2 = 1
        for c in cyc:
            cl2 *= 2 if c % 2 == 0 else 1
        I = cl2 + h
        pred = 0.5 * (I - 1) / (2 * h) * pi_sqrt
        obs = r["n_excl"] / (2 * h) - r["drawn_excl"]
        out.append({"P": P, "h": h, "Cl2": cl2, "observed_deficit": obs,
                    "predicted_mean_deficit": pred})
        obs_all.append(obs)
        pred_all.append(pred)
    obs_all, pred_all = np.array(obs_all), np.array(pred_all)
    summ = {
        "pi(sqrt N)": pi_sqrt,
        "P with positive deficit": int((obs_all > 0).sum()),
        "observed deficit: min, median, max": [float(obs_all.min()), float(np.median(obs_all)), float(obs_all.max())],
        "predicted mean deficit: min, median, max": [float(pred_all.min()), float(np.median(pred_all)), float(pred_all.max())],
        "sum observed": float(obs_all.sum()), "sum predicted": float(pred_all.sum()),
        "correlation": float(np.corrcoef(obs_all, pred_all)[0, 1]),
    }
    for o in out[:5] + [o for o in out if o["P"] in (11, 27, 58, 64, 89)]:
        print(o)
    print(summ)
    json.dump({"rows": out, "summary": summ}, open("bias.json", "w"), indent=1)


if __name__ == "__main__":
    main()
