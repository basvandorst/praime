"""Round 62 (preregister/round62.txt): information of the last k gap residues mod 30 about the next, data against the
hidden-residue model of round 61 at finite lam (q = 3, 5) and its lam -> infinity limit. usage: r62.py"""
import os, json, math, itertools
import numpy as np
import common as C
import r61

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
L0, L1 = 10 ** 10, 2 * 10 ** 10
KMAX = 4


def Hb(p):
    p = p[p > 0]; return float(-np.sum(p * np.log2(p)))


def mi_from_joint(J, k):
    """J: array of shape (15,)*(k+1), last axis = next gap; returns I(next; previous k)."""
    J = J / J.sum()
    past = J.sum(axis=-1); nxt = J.reshape(-1, 15).sum(axis=0)
    return Hb(past.ravel()) + Hb(nxt) - Hb(J.ravel())


def model(lam, k):
    c = np.arange(15)                            # g = 2c mod 30
    g = np.arange(2, 3000, 2)
    w = np.bincount((g // 2) % 15, weights=np.exp(-g / lam), minlength=15)
    n = k + 1
    J = np.zeros((15,) * n)
    for seq in itertools.product(range(15), repeat=n):
        pos = np.cumsum([0] + [2 * x for x in seq])
        n3 = 3 - len(set((pos % 3).tolist())); n5 = 5 - len(set((pos % 5).tolist()))
        if n3 <= 0 or n5 <= 0:
            continue
        J[seq] = np.prod(w[list(seq)]) * n3 * n5
    return mi_from_joint(J, k)


def data(k):
    J = np.zeros(15 ** (k + 1)); prev = None
    for s in range(L0, L1, 10 ** 8):
        P = np.nonzero(C.primes_mask(s, 10 ** 8))[0].astype(np.int64) + s
        if prev is not None:
            P = np.concatenate([prev, P])
        cc = ((np.diff(P) // 2) % 15).astype(np.int64)
        code = np.zeros(len(cc) - k, dtype=np.int64)
        for j in range(k + 1):
            code = code * 15 + cc[j:len(cc) - k + j]
        J += np.bincount(code, minlength=15 ** (k + 1))
        prev = P[-(k + 1):]
    N = J.sum(); cells = int(np.count_nonzero(J))
    Jr = J.reshape((15,) * (k + 1))
    I = mi_from_joint(Jr, k)
    # Miller-Madow on each entropy: I = H(past) + H(next) - H(joint); corrections (nz-1)/(2N ln2)
    nzp = int(np.count_nonzero(Jr.sum(axis=-1))); nzn = int(np.count_nonzero(Jr.reshape(-1, 15).sum(axis=0)))
    I += ((nzp - 1) + (nzn - 1) - (cells - 1)) / (2 * N * math.log(2))
    return I, int(N), cells


if __name__ == '__main__':
    lam = float(np.mean(np.log(np.linspace(L0, L1, 1001))))
    out = {'lambda': lam, 'k': {}}
    for k in range(1, KMAX + 1):
        loc = sum(r61.I_k_q(q, [k])[k] for q in (3, 5)) / math.log(2)
        mo = model(lam, k); dI, N, cells = data(k)
        out['k'][str(k)] = {'data': dI, 'model_finite_lambda': mo, 'limit_lambda_inf': loc, 'diff': dI - mo, 'n': N, 'cells': cells}
        print(k, {a: (round(b, 5) if isinstance(b, float) else b) for a, b in out['k'][str(k)].items()}, flush=True)
    out['saturation_bits'] = sum(((q - 2) / (q - 1)) * math.log((q - 1) / (q - 2)) for q in (3, 5)) / math.log(2)
    json.dump(out, open(os.path.join(OUT, 'r62.json'), 'w'), indent=1)
