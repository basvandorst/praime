"""Zeros of the L-functions of the curves of a registered set (degree-2 round), PARI lfunzeros(lfuncreate(E), TOP) at
the default subdivision (divz = 8), then the completeness check and repair. One .npy file per curve (positive
ordinates; the multiplicity of the central zero and the repair log go to ID.log.json), so a run can stop and continue: each call works until
its time budget is used and can be repeated.

Completeness (registered rule): s(t) = #{0 < gamma <= t} - theta_E(t)/pi + r/2 on 200 equally spaced heights in
[LO, HI] and at LO and HI (theta_E(t) = (t/2) ln N - t ln 2 pi + Im ln Gamma(1 + i t), r = number of zeros at 0), its
running mean over 21 consecutive heights (a missed close pair lowers s by 2 for the rest of the band; the raw s of a
degree-2 L-function reaches |s| = 1.5 by itself), and the largest unfolded gap in [LO, HI] (u = theta_E/pi).
Repair trigger: max|s| > 1.5, max|running mean| > 0.8 or a gap > 3.0. Then lfunzeros is run again with divz = 64
(twice at most, then divz = 256) on [t_c - 25, t_c + 6] around the first height t_c where |running mean| > 0.8, on
[t - 12, t + 6] around the first height where |s| > 1.5, and on [t_g - 6, t_g + 6] around the largest gap; zeros not
already present (|dt| > 1e-7) are added, the set is re-checked and everything is logged. Exclusion (deg2_test.py):
max|running mean| > 1.0 or max|s| > 2.0 after the repair; the exclusion is reported.

usage: deg2_zeros.py SET.json OUTDIR PART NPARTS [SECONDS]      compute (resumable)
       deg2_zeros.py SET.json OUTDIR check                       completeness table of all computed sets
Run as (ulimit -v 1500000; timeout 600 nice python -I deg2_zeros.py ...)."""
import sys, os, json, math, time
import numpy as np, cypari2
from scipy.special import loggamma

SMAX, GMAX = 1.5, 3.0


def theta_E(N, t):
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(2 * math.pi) + np.imag(loggamma(1 + 1j * t))


def check(g, N, r, lo, hi):
    """max|s|, max|running mean of s|, first height with |s| > 1.5, first centre with |running mean| > 0.8,
    largest unfolded gap in [lo, hi] and where, mean of s"""
    g = np.sort(g); t = np.concatenate([[lo], np.linspace(lo, hi, 200), [hi]])
    s = np.searchsorted(g, t, side='right') - theta_E(N, t) / math.pi + r / 2
    sm = np.convolve(s, np.ones(21) / 21, mode='valid'); tc = t[10:10 + len(sm)]
    b = g[(g >= lo) & (g <= hi)]; u = theta_E(N, b) / math.pi; d = np.diff(u); i = int(np.argmax(d))
    fs = float(t[np.argmax(np.abs(s) > 1.5)]) if (np.abs(s) > 1.5).any() else None
    fm = float(tc[np.argmax(np.abs(sm) > 0.8)]) if (np.abs(sm) > 0.8).any() else None
    return float(np.abs(s).max()), float(np.abs(sm).max()), fs, fm, float(d[i]), float(b[i]), float(np.mean(s))


def flagged(ms, msm, gap):
    return ms > SMAX or msm > 0.8 or gap > GMAX


def main():
    S = json.load(open(sys.argv[1])); out = sys.argv[2]; os.makedirs(out, exist_ok=True)
    lo, hi = S['band']; top = hi + 2.0
    curves = S['curves']
    lg = lambda cid: os.path.join(out, '%s.log.json' % cid)
    if sys.argv[3] == 'check':
        n = 0; nz = 0
        for c in curves:
            f = os.path.join(out, '%s.npy' % c['id'])
            if not os.path.exists(f): continue
            g = np.load(f); m = json.load(open(lg(c['id']))); ms, msm, fs, fm, gap, where, mean = check(g, c['N'], m['r'], lo, hi)
            excl = msm > 1.0 or ms > 2.0; n += 1; nz += int(((g >= lo) & (g <= hi)).sum())
            print('%-10s N=%5d r=%d zeros(band)=%5d max|s|=%.2f max|run.mean|=%.2f mean s=%+.2f gap=%.2f at %.2f added=%d %s' % (
                c['id'], c['N'], m['r'], int(((g >= lo) & (g <= hi)).sum()), ms, msm, mean, gap, where, m['added'],
                'EXCLUDED' if excl else ('flag (searched)' if flagged(ms, msm, gap) else '')))
        print('%d sets computed, %d zeros in the band' % (n, nz))
        return
    pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
    part, nparts = int(sys.argv[3]), int(sys.argv[4]); budget = float(sys.argv[5]) if len(sys.argv) > 5 else 540.0
    t0 = time.time(); done = 0
    for c in curves[part::nparts]:
        f = os.path.join(out, '%s.npy' % c['id'])
        if os.path.exists(f): continue
        if time.time() - t0 > budget: break
        t1 = time.time()
        L = pari.lfuncreate(pari.ellinit(c['ainv']))
        z = np.array([float(x) for x in pari.lfunzeros(L, top)])
        r = int((np.abs(z) < 1e-6).sum()); g = np.sort(z[z > 1e-6]); n0 = len(g); added = []; rounds = []
        for it in range(3):
            ms, msm, fs, fm, gap, where, mean = check(g, c['N'], r, lo, hi)
            if not flagged(ms, msm, gap): break
            divz = 64 if it < 2 else 256
            wins = []
            if fm is not None: wins.append((max(0.5, fm - 25), min(top, fm + 6)))
            if fs is not None: wins.append((max(0.5, fs - 12), min(top, fs + 6)))
            if gap > GMAX: wins.append((max(0.5, where - 6), min(top, where + 6)))
            if not wins: wins.append((lo, hi))
            new = []
            for a, b in wins:
                new += [float(x) for x in pari.lfunzeros(L, [a, b], divz)]
            new = [x for x in new if x > 1e-6 and np.min(np.abs(g - x)) > 1e-7]
            rounds.append({'divz': divz, 'windows': wins, 'added': new})
            if new:
                g = np.sort(np.concatenate([g, new])); added += new
        ms, msm, fs, fm, gap, where, mean = check(g, c['N'], r, lo, hi)
        np.save(f, g)
        m = {'N': c['N'], 'r': r, 'n0': n0, 'added': len(added), 'added_at': added, 'rounds': rounds,
                        'maxs': ms, 'max_runmean': msm, 'gap': gap, 'flag': bool(flagged(ms, msm, gap)),
                        'excluded': bool(msm > 1.0 or ms > 2.0), 'sec': round(time.time() - t1, 1)}
        json.dump(m, open(lg(c['id']), 'w')); done += 1
        print('%s N=%d r=%d zeros=%d added=%d max|s|=%.2f run.mean=%.2f gap=%.2f %.0fs' % (c['id'], c['N'], r, len(g), len(added), ms, msm, gap,
                                                                        time.time() - t1), flush=True)
    left = sum(not os.path.exists(os.path.join(out, '%s.npy' % c['id'])) for c in curves[part::nparts])
    print('part %d/%d: computed %d, left %d, %.0fs' % (part, nparts, done, left, time.time() - t0))


if __name__ == '__main__':
    main()
