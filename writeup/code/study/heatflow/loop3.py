"""Loop 3: flowed Haar matrices with a stable method (no polynomial root finder).
Coefficients by FFT of prod(1 + x z_i) on roots of unity; the flowed real function on a fine grid by FFT; zeros = sign changes."""
import sys, numpy as np
from scipy.stats import unitary_group
from loop2 import K_limit

def flowed_angles(z, t, over=256):
    N = len(z); M = 1
    while M < N + 1: M *= 2
    x = np.exp(2j * np.pi * np.arange(M) / M)
    P = np.prod(1 + x[:, None] * z[None, :], axis=1)
    e = np.fft.fft(P)[: N + 1] / M                       # e_j
    j = np.arange(N + 1)
    a = e * np.exp(-t * j * (N - j))
    G = over * M
    th = 2 * np.pi * np.arange(G) / G
    F = np.fft.ifft(a, n=G) * G                          # sum_j a_j e^{i j th}
    f = F * np.exp(-0.5j * N * th) / np.sqrt(e[N])
    im = np.max(np.abs(f.imag)) / np.max(np.abs(f.real))
    r = f.real
    s = np.nonzero(r * np.roll(r, -1) < 0)[0]
    r1 = r[(s + 1) % G]
    ang = (s + r[s] / (r[s] - r1)) * 2 * np.pi / G       # zeros of the flowed function; eigenangle = pi - angle, irrelevant here
    return np.sort(ang), im

if __name__ == "__main__":
    N = int(sys.argv[1]); n_mat = int(sys.argv[2])
    us = (0.3, 0.6, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0)
    ks = np.array([int(round(u * N)) for u in us])
    Us = unitary_group.rvs(N, size=n_mat, random_state=np.random.default_rng(5))
    zs = [np.linalg.eigvals(U) for U in Us]
    for tau in [float(a) for a in sys.argv[3:]]:
        vals = np.zeros((n_mat, len(ks))); sp = []; bad = 0; immax = 0
        for m, z in enumerate(zs):
            ang, im = flowed_angles(z, tau / N ** 2); immax = max(immax, im)
            if len(ang) != N: bad += 1; continue
            vals[m] = np.abs(np.exp(1j * np.outer(ks, ang)).sum(axis=1)) ** 2 / N
            d = np.diff(np.append(ang, ang[0] + 2 * np.pi)) * N / (2 * np.pi); sp.append(d)
        good = n_mat - bad; vals = vals[np.abs(vals).sum(axis=1) > 0]
        sp = np.array(sp)
        print(f"N = {N}, tau = {tau}: {good} of {n_mat} matrices with exactly N real zeros; imaginary part {immax:.1e}; "
              f"spacing variance {sp.var():.5f}, smallest spacing {sp.min():.4f}")
        print("    " + "  ".join(f"u={u:4.2f}: {v:7.4f}+-{e:6.4f}" + (f" ({sum(K_limit(u, tau)):6.4f})" if u < 1 else "")
                               for u, v, e in zip(us, vals.mean(axis=0), vals.std(axis=0) / np.sqrt(good))))
