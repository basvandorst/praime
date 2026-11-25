"""synth.py: synthetic |L(1/2+it)| on the data grid (t = 100 + k/20, k = 0..18000) for a discriminant D, in two
worlds (random model; used for the mechanism check and the dry run of freeze_test.py):
  field(t) = W_j(h + w*delta(t)) + sum_{p <= 31, p not | D} Y_p(t),   h = t - a_j  on interval j = [a_j, a_j + 2 pi)
  Y_p(t)    = -ln|1 - chi_D(p) p^(-1/2 - it)|              (exact Euler factor, coherent over all t)
  S_p(t)    = -(1/pi) arg(1 - chi_D(p) p^(-1/2 - it))      (its counting wave: N_present = N_missing + S_p)
  delta(t)  = (2 pi / L(t)) sum_p S_p(t),  L(t) = ln(|D| t / 2 pi)
  W_j       = CUE(N_j) log-modulus on the circle mapped onto interval j (N_j = round(L) at its centre), with the
              Fourier modes k < 4 removed (they are replaced by the explicit primes p <= 31, ln p < 3.5);
              independent per interval.
  world 'wave' (w = 1): the zeros carry the primes' counting waves (round 147's closure);
  world 'ghk'  (w = 0): the zeros part is independent of the Euler factors (GHK splitting).
usage: synth.py OUTDIR world seed D [D ...]"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modellib as ml
import chartools as ct

T0, H, NPT = 100.0, 0.05, 18001
THETA = 2 * np.pi
K0 = 4


def kron_any(D, p):
    return ct.kron(D, p)


def synth(D, world, rng):
    t = T0 + H * np.arange(NPT)
    q = abs(D)
    L = np.log(q * t / (2 * np.pi))
    Ysum = np.zeros(NPT); Ssum = np.zeros(NPT)
    for p in [int(x) for x in ml.primes_upto(31)]:
        c = kron_any(int(D), p)
        if c == 0:
            continue
        z = 1 - c * p ** -0.5 * np.exp(-1j * t * np.log(p))
        Ysum += -np.log(np.abs(z))
        Ssum += -np.angle(z) / np.pi
    delta = (2 * np.pi / L) * Ssum * (1.0 if world == 'wave' else 0.0)
    j = np.minimum(((t - T0) // THETA).astype(int), 142)
    W = np.zeros(NPT)
    for jj in range(143):
        sel = np.nonzero(j == jj)[0]
        a = T0 + jj * THETA
        N = int(round(np.log(q * (a + THETA / 2) / (2 * np.pi))))
        th = ml.cue_angles(rng, 1, N)[0]
        h = t[sel] - a + delta[sel]
        w = np.zeros(len(sel))
        for x in th:
            w += np.log(np.abs(2 * np.sin((h - x) / 2)) + 1e-300)
        for k in range(1, K0):
            Tk = np.exp(-1j * k * th).sum()
            w += np.real(np.exp(1j * k * h) * Tk) / k
        W[sel] = w
    return np.exp(W + Ysum)


def main():
    out, world, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
    os.makedirs(out, exist_ok=True)
    for D in [int(x) for x in sys.argv[4:]]:
        rng = np.random.default_rng([seed, abs(D), 1 if D > 0 else 0])
        np.save(os.path.join(out, 'absL_%d.npy' % D), synth(D, world, rng))


if __name__ == '__main__':
    main()
