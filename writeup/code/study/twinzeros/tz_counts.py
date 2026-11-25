"""tz_counts.py -- bounded, resumable, low-memory counts of primes and prime patterns on a fixed grid in ln x.

Grid: edges X_k = floor(exp(k/2048)) (integers, computed with mpmath at 40 digits, so the same on every machine).
A run covers   pairs  p in (PA, PB]   (twin-type pairs, triples, Sophie Germain pairs; weights are products of ln)
               primes p in (QA, QB]   (sum of ln p by residue class for the moduli QS; for the regressors)
Cell 0 of a run holds everything at or below the first edge of the run; cell j >= 1 holds (X_{k0+j-1}, X_{k0+j}].

Usage:
  python tz_counts.py RUNNAME PA PB QA QB [budget_seconds]
  discovery:  python tz_counts.py disc 0 2000000000 0 4000000000
  holdout:    python tz_counts.py hold 2000000000 20000000000 4000000000 40000000000
The run saves a state file RUNNAME_state.npz after each budget slice; rerun the same command until it prints DONE.
Output: RUNNAME.npz + RUNNAME.sha256. Arrays (per cell): th, cnt (sum of ln p, number of primes), cQ[:, r] (sum of ln p
over p = r mod Q), P{h}/N{h} (sum of ln p ln(p+h) / number of pairs with p, p+h prime), T026, T046 (triples p, p+2, p+6 and
p, p+4, p+6, weight ln p ln(p+a) ln(p+6)), SG (p, 2p+1, weight ln p ln(2p+1)); th_2, cQ_2: the same prime sums on the doubled
grid, cells (2 X_{j-1}, 2 X_j] (only for primes in (QA, QB]).
"""
import sys, os, time, hashlib
import numpy as np
import mpmath
from sieve_core import odd_flags, base_primes

DU_INV = 2048
QS = (3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 17, 19, 20, 21, 23, 24, 28, 29)
QS2 = (3, 4, 5, 7, 11, 12, 13, 15, 20, 21)   # moduli also kept on the doubled grid 2*X_k (for the Sophie Germain regressors)
HS = (2, 4, 6, 8, 10, 12, 30)
NSEG = 1 << 23          # odd numbers per segment (1.68e7 integers)
X_START = 10 ** 4       # first grid edge at or below 1e4


def edges_for(lo, hi):
    """Integer grid edges covering (lo, hi]: X_k for k0 <= k <= k1, X_k0 <= max(lo, X_START), X_k1 >= hi."""
    mpmath.mp.dps = 40
    base = max(lo, X_START)
    k0 = int(mpmath.floor(mpmath.log(base) * DU_INV))
    while int(mpmath.floor(mpmath.exp(mpmath.mpf(k0) / DU_INV))) > base:
        k0 -= 1
    k1 = int(mpmath.ceil(mpmath.log(hi) * DU_INV))
    while int(mpmath.floor(mpmath.exp(mpmath.mpf(k1) / DU_INV))) < hi:
        k1 += 1
    ks = np.arange(k0, k1 + 1)
    X = np.array([int(mpmath.floor(mpmath.exp(mpmath.mpf(int(k)) / DU_INV))) for k in ks], dtype=np.int64)
    return ks, X


def main():
    run = sys.argv[1]
    PA, PB, QA, QB = (int(float(a)) for a in sys.argv[2:6])
    budget = float(sys.argv[6]) if len(sys.argv) > 6 else 520.0
    t_start = time.time()
    state = run + "_state.npz"
    lo, hi = min(PA, QA), max(PB, QB)
    ks, X = edges_for(lo, hi)
    nc = len(X)                       # cells 0..nc-1 ; cell j>=1 is (X[j-1], X[j]]
    if os.path.exists(state):
        z = np.load(state)
        A = {k: z[k].copy() for k in z.files if k not in ("next_s", "ks", "X", "meta")}
        s = int(z["next_s"])
        assert np.array_equal(z["X"], X) and list(z["meta"]) == [PA, PB, QA, QB]
        print("resume at", s)
    else:
        A = {"th": np.zeros(nc), "cnt": np.zeros(nc), "th_2": np.zeros(nc)}
        for q in QS:
            A["c%d" % q] = np.zeros((nc, q))
        for q in QS2:
            A["c%d_2" % q] = np.zeros((nc, q))
        for h in HS:
            A["P%d" % h] = np.zeros(nc)
            A["N%d" % h] = np.zeros(nc)
        for t in ("T026", "T046", "SG"):
            A[t] = np.zeros(nc)
            A["N" + t] = np.zeros(nc)
        s = (lo // 2) * 2 + 1         # first odd number > lo-1
        # the prime 2
        if lo < 2 <= hi:
            j = int(np.searchsorted(X, 2, side="left"))
            if QA < 2 <= QB:
                A["th"][j] += np.log(2.0); A["cnt"][j] += 1
                for q in QS:
                    A["c%d" % q][j, 2 % q] += np.log(2.0)
                j2 = int(np.searchsorted(2 * X, 2, side="left"))
                A["th_2"][j2] += np.log(2.0)
                for q in QS2:
                    A["c%d_2" % q][j2, 2 % q] += np.log(2.0)
            if PA < 2 <= PB:          # (2, 5) is a Sophie Germain pair; no even-gap pair or triple starts at 2
                A["SG"][j] += np.log(2.0) * np.log(5.0); A["NSG"][j] += 1
    sp = base_primes(int((2 * hi + 64) ** 0.5) + 2)[1:]
    t_last_save = time.time()
    m840 = {q: (np.arange(840) % q).astype(np.int64) for q in QS if 840 % q == 0}
    while s <= hi:
        n = NSEG
        f = odd_flags(s, n + 16, sp)
        idx = np.flatnonzero(f[:n])
        p = s + 2 * idx.astype(np.int64)
        keep = p <= hi
        idx, p = idx[keep], p[keep]
        lnp = np.log(p.astype(np.float64))
        j = np.searchsorted(X, p, side="left")
        assert j.max(initial=0) < nc
        # primes: class sums
        mq = (p > QA) & (p <= QB)
        if mq.any():
            pq, jq, lq = p[mq], j[mq], lnp[mq]
            jmin, jmax = int(jq.min()), int(jq.max())
            span = jmax - jmin + 1
            jj = jq - jmin
            A["th"][jmin:jmax + 1] += np.bincount(jj, weights=lq, minlength=span)
            A["cnt"][jmin:jmax + 1] += np.bincount(jj, minlength=span)
            r840 = pq % 840
            for q in QS:
                r = m840[q][r840] if q in m840 else pq % q
                A["c%d" % q][jmin:jmax + 1] += np.bincount(jj * q + r, weights=lq, minlength=span * q).reshape(span, q)
            # the doubled grid: cells (2 X_{j-1}, 2 X_j]; only primes up to 2 X[-1] fit
            j2 = np.searchsorted(2 * X, pq, side="left")
            ok = j2 < nc
            if ok.any():
                j2, pq2, lq2 = j2[ok], pq[ok], lq[ok]
                jmin2, jmax2 = int(j2.min()), int(j2.max())
                span2 = jmax2 - jmin2 + 1
                jj2 = j2 - jmin2
                A["th_2"][jmin2:jmax2 + 1] += np.bincount(jj2, weights=lq2, minlength=span2)
                r840b = pq2 % 840
                for q in QS2:
                    r = m840[q][r840b] if q in m840 else pq2 % q
                    A["c%d_2" % q][jmin2:jmax2 + 1] += np.bincount(jj2 * q + r, weights=lq2, minlength=span2 * q).reshape(span2, q)
        # pairs, triples, Sophie Germain
        mp_ = (p > PA) & (p <= PB)
        if mp_.any():
            ip, pp, jp, lp = idx[mp_], p[mp_], j[mp_], lnp[mp_]
            jmin, jmax = int(jp.min()), int(jp.max())
            span = jmax - jmin + 1
            jj = jp - jmin
            for h in HS:
                g = f[ip + h // 2]
                w = lp[g] * np.log((pp[g] + h).astype(np.float64))
                A["P%d" % h][jmin:jmax + 1] += np.bincount(jj[g], weights=w, minlength=span)
                A["N%d" % h][jmin:jmax + 1] += np.bincount(jj[g], minlength=span)
            for name, a, b in (("T026", 1, 3), ("T046", 2, 3)):
                g = f[ip + a] & f[ip + b]
                w = lp[g] * np.log((pp[g] + 2 * a).astype(np.float64)) * np.log((pp[g] + 2 * b).astype(np.float64))
                A[name][jmin:jmax + 1] += np.bincount(jj[g], weights=w, minlength=span)
                A["N" + name][jmin:jmax + 1] += np.bincount(jj[g], minlength=span)
            fs = odd_flags(2 * s + 1, 2 * n, sp)          # odd numbers 2s+1, 2s+3, ...; 2p+1 sits at index 2*idx
            g = fs[2 * ip]
            w = lp[g] * np.log((2 * pp[g] + 1).astype(np.float64))
            A["SG"][jmin:jmax + 1] += np.bincount(jj[g], weights=w, minlength=span)
            A["NSG"][jmin:jmax + 1] += np.bincount(jj[g], minlength=span)
            del fs
        s += 2 * n
        if time.time() - t_start > budget and s <= hi:
            np.savez(state + ".tmp.npz", next_s=s, ks=ks, X=X, meta=np.array([PA, PB, QA, QB]), **A)
            os.replace(state + ".tmp.npz", state)
            print("paused at s = %d (%.1f%% done); rerun to continue" % (s, 100 * (s - lo) / (hi - lo)))
            return 3
    out = run + ".npz"
    np.savez(out, ks=ks, X=X, meta=np.array([PA, PB, QA, QB]), **A)
    h = hashlib.sha256(open(out, "rb").read()).hexdigest()
    open(run + ".sha256", "w").write(h + "  " + out + "\n")
    if os.path.exists(state):
        os.remove(state)
    print("DONE", out, h, "%.1fs" % (time.time() - t_start))
    return 0


if __name__ == "__main__":
    sys.exit(main())
