"""Primes from bands of zeta zeros: Landau's formula with an exact window decoder.

A band of zeros gamma in [c - W/2, c + W/2], weighted by a Tukey window w (flat, with cosine
fades of length W/7 at both ends), gives after demodulation by e^{-i c ln x}

  r(x) = -2 pi sqrt(x) / What(0) * sum_gamma w(gamma) e^{i (gamma - c) ln x}
       = sum_m Lambda(m) sqrt(x/m) e^{-i c ln m} What(ln(x/m)) / What(0),

where What is the Fourier transform of w (real and even) and m runs over the prime powers. So a
band is a lowpass view of the primes, at resolution about 2 pi x / W numbers, in which every
prime power m carries the known phase e^{-i c ln m}. The lowpass [0, T] is the special case c = 0
with the mirror zeros: r(x) = 2 pi sqrt(x) (sqrt(x) - 2 sum w cos(gamma ln x)) / What(0).

Decoder (as transformation 134 of prime.html, on complex samples; with --M the M-algorithm of
RIEMANN.md section 8: the M best decision histories by summed window cost, a decision committed
when it is 24 candidates old): walk up over the candidates
(odd numbers, or numbers coprime to the wheel); the echo of every decided prime power is
subtracted from all samples once; the next k candidates are searched exhaustively (0 or ln m);
the following candidates (12 envelope scales) are coins with the mean prime density (1/rho)/ln m,
whose covariance whitens the samples; keep the first decision and slide. Samples are taken at a third of the
envelope scale 2 pi n / W, from 4.4 to 3 scales around the window, with noise 0.01 per sample
scaled by the sampling step.

With --cold S nothing is known below S - 6 scales: the numbers there are coins as well.

usage:
  python3 band_decoder.py SPEC N [--k 12] [--M 8] [--wheel 2,3,5,7] [--cold S [--warm]] [--out file]
  SPEC: '+'-joined list of
        Odlyzko table (writeup/data/riemann/odlyzko_zeros3.txt, ...),
        low:T:zeros.f64        zeros up to T (lowpass),
        band:T1:W:zeros.f64    zeros in (T1, T1 + W].
Prints the first wrong decisions (checked against a sieve) and writes the decided primes.
"""
import argparse
import math
import re
import time

import mpmath as mp
import numpy as np

mp.mp.dps = 50
FADE = 1 / 7
OVERS, SPAN_L, SPAN_R, SIG, FAR = 3.0, 4.4, 3.0, 0.01, 12.0


def tukey_ft(v, A, L):
    """Fourier transform of the Tukey window of width A with cosine fades of length L."""
    v = np.asarray(v, float)
    a = A - L
    x = a * v / 2
    s = np.where(np.abs(x) < 1e-12, 1.0, np.sin(x) / np.where(x == 0, 1, x))
    den = 1 - (v * L / np.pi) ** 2
    small = np.abs(den) < 1e-9
    h = np.where(small, np.pi / 4, np.cos(v * L / 2) / np.where(small, 1, den))
    return a * s * h


class Band:
    def __init__(self, d, base=0, low=False):
        d = np.sort(np.asarray(d, float))
        self.low = low
        if low:
            T = d[-1]
            self.A, self.L = 2 * T, T * FADE
            w = np.ones_like(d)
            t = d > T - self.L
            w[t] = 0.5 * (1 + np.cos(np.pi * (d[t] - (T - self.L)) / self.L))
            self.g, self.w, self.c = d, w, mp.mpf(0)
        else:
            lo, hi = d[0], d[-1]
            self.A, self.L = hi - lo, (hi - lo) * FADE
            w = np.ones_like(d)
            t = d > hi - self.L
            w[t] = 0.5 * (1 + np.cos(np.pi * (d[t] - (hi - self.L)) / self.L))
            t = d < lo + self.L
            w[t] = 0.5 * (1 - np.cos(np.pi * (d[t] - lo) / self.L))
            dc = (lo + hi) / 2
            self.g, self.w = d - dc, w
            self.c = mp.mpf(base) + mp.mpf(dc)
        self.W0 = float(tukey_ft(0.0, self.A, self.L))
        self.width = self.A / 2 if low else self.A
        self.count = len(d)

    def samples(self, xs):
        xs = np.asarray(xs, float)
        u = np.log(xs)
        out = np.zeros(len(xs), complex)
        for a in range(0, len(xs), 256):
            ph = np.outer(u[a:a + 256], self.g)
            if self.low:
                S = 2 * (np.cos(ph) @ self.w)
                out[a:a + 256] = (np.sqrt(xs[a:a + 256]) - S) * 2 * np.pi
            else:
                out[a:a + 256] = -2 * np.pi * ((np.cos(ph) @ self.w) + 1j * (np.sin(ph) @ self.w))
        return out * np.sqrt(xs) / self.W0

    def phases(self, M):
        """e^{-i c ln m} for m = 0..M (0 at m = 0, 1)."""
        out = np.zeros(M + 1, complex)
        if self.low:
            out[1:] = 1
            return out
        two = 2 * mp.pi
        for m in range(2, M + 1):
            out[m] = complex(mp.expj(-mp.fmod(self.c * mp.log(m), two)))
        return out

    def kern(self, xs, ms):
        X = np.asarray(xs, float)[:, None]
        M = np.asarray(ms, float)[None, :]
        return np.sqrt(X / M) * tukey_ft(np.log(X / M), self.A, self.L) / self.W0

    def scale(self, n):
        """envelope scale in numbers at n: 2 pi n / width."""
        return 2 * np.pi * n / self.width if not self.low else np.pi * n / self.width


def read_odlyzko(path):
    txt = open(path).read().splitlines()
    base = int(re.search(r'gamma - (\d+)', ' '.join(txt[:6])).group(1))
    d = np.array([float(s) for s in txt if re.fullmatch(r'\s*\d+\.\d+\s*', s)])
    return Band(d, base)


def parse_spec(spec):
    bands, cache = [], {}
    for part in spec.split('+'):
        f = part.split(':')
        if f[0] == 'low':
            G = cache.setdefault(f[2], np.fromfile(f[2]))
            bands.append(Band(G[G <= float(f[1])], low=True))
        elif f[0] == 'band':
            G = cache.setdefault(f[3], np.fromfile(f[3]))
            T1, W = float(f[1]), float(f[2])
            bands.append(Band(G[(G > T1) & (G <= T1 + W)]))
        else:
            bands.append(read_odlyzko(part))
    return bands


def sieve(N):
    s = np.ones(N + 1, bool)
    s[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return s


def decode(bands, N, k=12, wheel=(2,), cold=None, stop_after=3, burn=6.0, pastr=0.12, warm=False, M=1, D=24):
    Xmax = N + 4000
    isp = sieve(Xmax)
    lam = np.zeros(Xmax + 1)
    root = np.zeros(Xmax + 1, np.int64)
    for p in np.nonzero(isp)[0]:
        q = int(p)
        while q <= Xmax:
            lam[q], root[q] = math.log(p), p
            q *= int(p)
    perfect = np.zeros(Xmax + 1, bool)
    for e in range(2, 20):
        r = 2
        while r ** e <= Xmax:
            perfect[r ** e] = True
            r += 1
    rho = float(np.prod([1 - 1 / q for q in wheel]))
    cand = np.ones(Xmax + 1, bool)
    cand[:2] = False
    for q in wheel:
        cand[::q] = False
    cand &= ~perfect
    xs = np.arange(2, Xmax + 1)
    res = [b.samples(xs) for b in bands]
    ph = [b.phases(Xmax) for b in bands]

    def subtract(m, value):
        for bi, b in enumerate(bands):
            res[bi] -= value * ph[bi][m] * b.kern(xs, [m])[:, 0]

    start = 3
    if cold:
        start = max(3, int(cold - burn * max(b.scale(cold) for b in bands)))
    # known: prime powers of the wheel primes
    for m in range(2, Xmax + 1):
        if lam[m] > 0 and root[m] in wheel:
            subtract(m, lam[m])
    later_pp, pp_coins = {}, []
    for m in range(start, Xmax + 1):
        if lam[m] > 0 and perfect[m] and root[m] not in wheel:
            if root[m] >= start:
                later_pp.setdefault(int(root[m]), []).append(m)
            else:
                pp_coins.append(m)
    S0 = start
    if cold and warm:   # oracle: the true values between S0 and S are given
        for m in range(S0, cold):
            if root[m] == m and m not in wheel:
                subtract(m, lam[m])
                for q in later_pp.pop(m, []):
                    subtract(q, lam[q])
        start = cold
    allc = [m for m in range(start, Xmax + 1) if cand[m]]
    past = np.array([m for m in range(max(3, int(S0 * (1 - pastr))), S0)
                     if cand[m] or (perfect[m] and lam[m] > 0 and root[m] not in wheel)], int) if cold else np.array([], int)
    pats = None
    errors, decided = [], []
    hist = [([], 0.0)]   # M-algorithm: (uncommitted decisions, summed window cost)
    stopped = False

    def commit(mm, dec):
        if dec:
            subtract(mm, math.log(mm))
            for q in later_pp.get(mm, []):
                subtract(q, math.log(mm))
        if cold and mm < cold:
            return False
        if dec:
            decided.append(mm)
        if dec != bool(isp[mm]):
            errors.append(mm)
            return bool(stop_after and len(errors) >= stop_after)
        return False

    for i, n in enumerate(allc):
        if n > N or stopped:
            break
        U = np.array(allc[i:i + k])
        far_n = int(max(60, FAR * max(b.scale(n) for b in bands) * rho))
        F = np.array(allc[i + k:i + k + far_n], int)
        if cold:
            near = [m for m in pp_coins if abs(m - n) < 4 * max(b.scale(n) for b in bands) + 50]
            F = np.concatenate([past, F, np.array(near, int)])
        pF = np.minimum(0.95, (1 / rho) / np.log(F))
        vF = pF * (1 - pF) * np.log(F) ** 2
        ys, As, Ks, sg, SX = [], [], [], [], []
        for bi, b in enumerate(bands):
            sc = b.scale(n)
            step = max(1, int(sc / OVERS))
            sx = np.arange(max(2, n - int(SPAN_L * sc)), U[-1] + int(SPAN_R * sc) + 5, step)
            SX.append(sx)
            r = res[bi][sx - 2]
            Kf = b.kern(sx, F) * ph[bi][F][None, :]
            r = r - Kf @ (pF * np.log(F))
            Ku = b.kern(sx, U) * (np.log(U) * ph[bi][U])[None, :]
            if b.low:
                ys.append(r.real); As.append(Ku.real); Ks.append(Kf.real)
            else:
                ys.append(np.concatenate([r.real, r.imag]))
                As.append(np.vstack([Ku.real, Ku.imag]))
                Ks.append(np.vstack([Kf.real, Kf.imag]))
            sg.append(np.full(len(ys[-1]), SIG / step))
        y, A, K = np.concatenate(ys), np.vstack(As), np.vstack(Ks)
        C = np.diag(np.concatenate(sg)) + (K * vF[None, :]) @ K.T
        Lc = np.linalg.cholesky(C)
        Q, R = np.linalg.qr(np.linalg.solve(Lc, A))
        if pats is None or pats.shape[1] != len(U):
            pats = ((np.arange(2 ** len(U))[:, None] >> np.arange(len(U))[None, :]) & 1).astype(float)
        if M == 1:
            z = Q.T @ np.linalg.solve(Lc, y)
            best = pats[np.argmin(((z[None, :] - pats @ R.T) ** 2).sum(1))]
            stopped = commit(n, bool(best[0]))
            continue
        # M-algorithm: whitening and QR are shared; each history differs only by its uncommitted echo
        P = pats @ R.T
        children = []
        for unc, met in hist:
            yh = y
            ones = [m for m, bit in unc if bit]
            if ones:
                adj = []
                for bi, b in enumerate(bands):
                    col = (b.kern(SX[bi], ones) * (np.log(ones) * ph[bi][ones])[None, :]).sum(1)
                    adj.append(col.real if b.low else np.concatenate([col.real, col.imag]))
                yh = y - np.concatenate(adj)
            cost = ((Q.T @ np.linalg.solve(Lc, yh))[None, :] - P) ** 2
            cost = cost.sum(1)
            children.append((unc + [(n, False)], met + cost[pats[:, 0] == 0].min()))
            children.append((unc + [(n, True)], met + cost[pats[:, 0] == 1].min()))
        children.sort(key=lambda h: h[1])
        hist = children[:M]
        if len(hist[0][0]) > D:
            m0, b0 = hist[0][0][0]
            hist = [(u[1:], met) for u, met in hist if u[0] == (m0, b0)]
            stopped = commit(m0, b0)
    if M > 1 and not stopped:
        for m0, b0 in hist[0][0]:
            if commit(m0, b0):
                break
    return decided, errors


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('spec')
    ap.add_argument('N', type=int)
    ap.add_argument('--k', type=int, default=12)
    ap.add_argument('--wheel', default='2')
    ap.add_argument('--cold', type=int, default=None)
    ap.add_argument('--stop', type=int, default=3)
    ap.add_argument('--out', default=None)
    ap.add_argument('--warm', action='store_true', help='oracle: true values between S0 and S given')
    ap.add_argument('--M', type=int, default=1, help='decision histories kept (M-algorithm)')
    a = ap.parse_args()
    t0 = time.time()
    bands = parse_spec(a.spec)
    wheel = tuple(int(q) for q in a.wheel.split(','))
    decided, errors = decode(bands, a.N, k=a.k, wheel=wheel, cold=a.cold, stop_after=a.stop, warm=a.warm, M=a.M)
    width = sum(b.width for b in bands)
    nz = sum(b.count for b in bands)
    reach = errors[0] if errors else a.N
    centres = ', '.join(f'{float(b.c):.6g}' for b in bands)
    print(f"{len(bands)} band(s) at centre {centres}, total width {width:.0f}, {nz} zeros, wheel {wheel}"
          + (f", cold start at {a.cold}" if a.cold else "") + (f", M = {a.M}" if a.M > 1 else "")
          + f": {len(decided)} primes decided; first wrong {errors[:a.stop]}; "
          f"{'reach' if errors else 'clean to'} {reach} = {reach / width:.2f} x width [{time.time() - t0:.0f} s]")
    if a.out:
        with open(a.out, 'w') as f:
            f.write(f"# {a.spec} N={a.N} wheel={a.wheel} cold={a.cold}; first wrong {errors[:a.stop]}\n")
            f.write('\n'.join(str(p) for p in decided) + '\n')
