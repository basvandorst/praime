"""Growth of LHS / [(MT + N min(M,T)) sum|c|^2] (Stucky Lemma 4.2, T <= M) with N, dyadic support (N/2, N].

(a) c = 1 on rational ideals (k), k^2 in (N/2, N]   -> predicted growth ~ sqrt(N)  (the counterexample)
(b) coefficients of the type that occur: 1, mu, d, and the convolutions 1*mu, 1*1 of two dyadic pieces
    -> the content weight sum |c|^2 cont / sum |c|^2 stays O(log N), so the ratio should grow at most like log N.
T = 1, M = 10^5 (M >> N), the regime in which the form fails.
"""
import sys
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ideals import enumerate_ideals
from mv_test import kernel, lhs


def normalize(a, b):
    """first-quadrant generator (a >= 1, b >= 0) of the ideal (a + bi)"""
    for _ in range(4):
        if a > 0 and b >= 0:
            return a, b
        a, b = -b, a  # multiply by i
    raise ValueError


def convolution(Nout, N1, N2, f2):
    """c(a) = sum_{b c = a, Nb in (N1/2, N1], Nc in (N2/2, N2]} f2(c) as dict on generators"""
    I1 = enumerate_ideals(N1, N1 // 2)
    I2 = enumerate_ideals(N2, N2 // 2)
    out = {}
    v2 = f2(I2)
    for i in range(len(I1['a'])):
        a1, b1 = int(I1['a'][i]), int(I1['b'][i])
        for j in range(len(I2['a'])):
            if v2[j] == 0:
                continue
            a2, b2 = int(I2['a'][j]), int(I2['b'][j])
            g = normalize(a1 * a2 - b1 * b2, a1 * b2 + a2 * b1)
            out[g] = out.get(g, 0) + v2[j]
    return out


def main():
    M, T = 100000, 1.0
    print("(a) rational ideals (k), k^2 in (N/2, N], T = 1, M = 10 N: ratio / sqrt(N)")
    for N in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
        k = np.array([k for k in range(1, math.isqrt(N) + 1) if N // 2 < k * k <= N], dtype=float)
        u = 2 * (np.log(k)[:, None] - np.log(k)[None, :])
        KT = 2 * T * np.sinc(T * u / math.pi)
        Ma = 10 * N  # the form fails when M >> sqrt(N); take M = 10 N
        L = (2 * Ma + 1) * KT.sum()
        r = L / ((Ma * T + N * min(Ma, T)) * len(k))
        print(f"   N = {N:>9}: #c = {len(k):5d}  ratio = {r:9.3f}   ratio/sqrt(N) = {r / math.sqrt(N):.4f}")
    print("\n(b) coefficient types of the method, dyadic support (N/2, N]")
    for N in (1000, 2000, 4000, 6000):
        I = enumerate_ideals(N, N // 2)
        idx = np.arange(len(I['a']))
        K = kernel(I, idx, M, T)
        pos = {(int(I['a'][i]), int(I['b'][i])): i for i in idx}
        fams = {'rational (k)': (I['b'] == 0).astype(float), 'ones': np.ones(len(idx)),
                'Moebius': I['mu'].astype(float), 'd(a)': I['d'].astype(float)}
        # convolutions of two dyadic pieces with N1 N2 ~ N
        N1 = int(round(math.sqrt(N))) * 2
        N2 = N // N1 * 2
        for name, f2 in (('1*mu', lambda J: J['mu'].astype(float)), ('1*1', lambda J: np.ones(len(J['a'])))):
            conv = convolution(N, N1, N2, f2)
            c = np.zeros(len(idx))
            for g, v in conv.items():
                if g in pos:
                    c[pos[g]] = v
            fams[f'{name} (N1={N1},N2={N2})'] = c
        line = []
        for name, c in fams.items():
            c = c.astype(complex)
            s2 = float(np.sum(np.abs(c) ** 2))
            if s2 == 0:
                continue
            L = lhs(K, c)
            cw = float(np.sum(np.abs(c) ** 2 * I['cont'])) / s2
            line.append(f"{name}: {L / ((M * T + N * min(M, T)) * s2):.3f} (cont-weight {cw:.2f})")
        print(f"   N = {N}: " + ";  ".join(line))
        print(f"            log N = {math.log(N):.2f}, sqrt(N) = {math.sqrt(N):.1f}")


if __name__ == '__main__':
    main()
