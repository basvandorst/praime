"""J_inf = lim I(g_{n+1}; g_{n-1} | g_n): sum over odd primes of the local averages over (a, b, c) mod p of
f4 ln(f4 f2(b) / (f3(a,b) f3(b,c))), f_k the local factor (1 - nu/p)/(1 - 1/p)^k of the residue set; p = 2 gives 0.
Exact enumeration for p <= 400, a fitted c/p^4 + d/p^5 tail beyond (checked on 200 < p <= 400)."""
import math, json, os
import numpy as np
import common as C


def local(p):
    b, c = np.meshgrid(np.arange(p), np.arange(p), indexing='ij'); b = b.ravel(); c = c.ravel()
    f = lambda n, k: (1 - n / p) / (1 - 1 / p) ** k
    def nu(*pts):
        R = np.stack([x % p for x in pts], 1); R.sort(1)
        return 1 + np.count_nonzero(np.diff(R, axis=1), axis=1)
    z = np.zeros_like(b)
    f2b = f(nu(z, b), 2); f3bc = f(nu(z, b, b + c), 3)
    tot = 0.0
    for a in range(p):
        av = np.full_like(b, a)
        f4 = f(nu(z, av, av + b, av + b + c), 4); f3ab = f(nu(z, av, av + b), 3)
        m = f4 > 0
        tot += float(np.sum(f4[m] * np.log(f4[m] * f2b[m] / (f3ab[m] * f3bc[m]))))
    return tot / p ** 3


if __name__ == '__main__':
    ps = [int(q) for q in C.small_primes(400)[1:]]
    J = {p: local(p) for p in ps}
    head = sum(J.values())
    big = np.array([p for p in ps if p > 200], float); vals = np.array([J[int(p)] for p in big])
    A = np.stack([big ** -4.0, big ** -5.0], 1); coef = np.linalg.lstsq(A, vals, rcond=None)[0]
    P = C.small_primes(20_000_000).astype(float); P = P[P > 400]
    tail = float(np.sum(coef[0] * P ** -4 + coef[1] * P ** -5))
    Jinf = head + tail
    out = {'J_3': J[3], 'J_5': J[5], 'J_7': J[7], 'head_p_le_400': head, 'tail': tail, 'J_inf_nats': Jinf, 'J_inf_bits': Jinf / math.log(2),
           'fit_check_max_rel': float(np.max(np.abs(A @ coef - vals) / np.abs(vals)))}
    print({k: (round(v, 8) if isinstance(v, float) else v) for k, v in out.items()})
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study', 'r58_theory.json'), 'w'), indent=1)
