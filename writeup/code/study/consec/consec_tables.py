"""consec_tables.py: transition tables of the labels of consecutive primes, per segment (bounded, resumable, low memory).

usage: consec_tables.py OUTDIR X0 X1 [--seg 5e6] [--cases C23,C31,C47,C39,G,E37] [--max-seconds 540]
  Primes p_n with X0 < p_n < X1 (X0, X1 multiples of the segment length).  A pair (p_n, p_(n+k)) is counted in the
  segment of p_n and only if p_(n+k) < X1 (nothing beyond X1 is ever computed).  One file OUTDIR/seg_<L>.npz per segment;
  existing files are skipped (resumable).  Stops cleanly after --max-seconds.

Per case ("schemes"; see consec_lib.py for the label codes):
  lab        label code per prime (alphabet A), alpha = abelian projection of the label
  mask       primes that carry the non-abelian functions f (the "fibre"), fib = class code inside the fibre (C classes)
Stored per segment and scheme:
  X{k}       A x A counts of (lab_n, lab_(n+k)), k = 1, 2, 3;  X3t: A x A x A triples (n, n+1, n+2)
  Xg         lag-1 table by gap index min(g/2, 100): (101, A, A)
  H{k}_*     mask-mask pair sums at lag k (consecutive in the full prime sequence): n, sf (F), sg (F), sfg (F x F),
             sff (F), sgg (F);   M{k}_*: the same for lag k inside the mask subsequence (next fibre prime)
  G1_*       lag-1 mask-mask sums by gap index: n, sfg (diagonal, F), sf, sg, sff, sgg
  B1_*       sums of f(p_n) over mask p_n by alpha(p_(n+k)) and of f(p_(n+k)) over mask p_(n+k) by alpha(p_n), k = 1, 2, 3
  runs       run-length histogram (C x 41) of the fibre label in the mask subsequence (edge runs dropped)
  D{H}_*     block dispersion, blocks of length H = 1e2 .. 1e6: per function f: number of blocks, sum T^2/V, sum (T^2/V)^2,
             sum T^2, sum V, sum T^4, sum V^2, sum T^2 V, number of blocks incl. empty ones, with T = sum_block f,
             V = sum_block f^2; the last row is the abelian sign chi (+1 / -1; zero rows for E37s)
  W4         per block of 1e4: C x C table of the fibre labels of consecutive mask-mask pairs (n, n+1)
  W4a        per block of 1e4: 2 x 2 table (abelian sign chi of p_n, p_(n+1)) (dihedral and G only)
"""
import sys, os, time, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from consec_lib import segment, F_DIH

GCAP = 100; HS = (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6); EXT = 4000; RUNCAP = 40
S2 = math.sqrt(2.0)

def schemes(seg):
    """dict name -> (lab int8, A, alpha (int array over label codes), mask bool, fib int (class in fibre), C, F list of float arrays, chi)"""
    p = seg['p']; out = {}
    for c in ('C23', 'C31', 'C47', 'C39'):
        if c not in seg: continue
        D = int(c[1:]); lab = seg[c]; A = 4 if D in (47, 39) else 3
        if D == 39:
            alpha = np.array([0, 1, 1, 2]); mask = (lab == 1) | (lab == 2); fib = np.where(lab == 2, 1, 0)
            F = [np.where(lab == 1, 1.0, np.where(lab == 2, -1.0, 0.0))]; C = 2
        else:
            alpha = np.array([0] + [1] * (A - 1)); mask = lab >= 1; fib = np.where(mask, lab - 1, 0); C = A - 1
            F = [np.where(mask, fv[np.clip(lab, 0, None)], 0.0) for fv in F_DIH[D]]
        chi = np.where(lab < 0, 0, np.where(lab == 0, -1, 1)).astype(np.float64)
        out[c] = (lab, A, alpha, mask, fib, C, F, chi)
    if 'G' in seg:
        lab = seg['G']; th = np.nan_to_num(seg['G_th']); mask = lab >= 1
        F = [np.where(mask, S2 * np.cos(4 * m * th), 0.0) for m in (1, 2)]
        chi = np.where(lab < 0, 0, np.where(lab == 0, -1, 1)).astype(np.float64)
        out['G'] = (lab, 5, np.array([0, 1, 1, 1, 1]), mask, np.where(mask, lab - 1, 0), 4, F, chi)
    if 'E37_ap' in seg:
        ap = seg['E37_ap']; ok = (p != 37)
        u1 = ap / np.sqrt(p.astype(np.float64))
        # sign of a_p
        lab = np.where(~ok, -1, np.where(ap > 0, 1, np.where(ap < 0, 0, 2))).astype(np.int8); mask = ok & (ap != 0)
        F = [np.where(mask, np.sign(ap).astype(np.float64), 0.0), np.where(mask, u1, 0.0), np.where(mask, u1 * u1 - 1, 0.0)]
        out['E37s'] = (lab, 3, np.array([0, 0, 0]), mask, np.where(mask, lab, 0), 2, F, np.zeros(len(p)))
        # a_p mod 2 with the abelian part (p/37)
        leg = np.zeros(37, np.int8); leg[(np.arange(1, 37) ** 2) % 37] = 1; leg = np.where(leg == 1, 1, -1); leg[0] = 0
        l37 = leg[p % 37]
        lab = np.where(~ok, -1, np.where(l37 == -1, 0, np.where(ap % 2 == 0, 1, 2))).astype(np.int8); mask = lab >= 1
        F = [np.where(lab == 1, S2, np.where(lab == 2, -1 / S2, 0.0))]
        chi = np.where(lab < 0, 0, np.where(lab == 0, -1, 1)).astype(np.float64)
        out['E37m2'] = (lab, 3, np.array([0, 1, 1]), mask, np.where(mask, lab - 1, 0), 2, F, chi)
        # a_p mod 3 with the abelian part p mod 3
        ok3 = ok & (p != 3); r3 = (p % 3 == 2).astype(np.int64); t = ap % 3
        lab = np.where(~ok3, -1, 3 * r3 + t).astype(np.int8); mask = ok3
        P0 = np.where(r3 == 1, 0.5, 0.25); s1 = np.sqrt(P0 * (1 - P0)); s2 = np.where(r3 == 1, math.sqrt(0.5), math.sqrt(0.75))
        F = [np.where(mask, ((t == 0) - P0) / s1, 0.0), np.where(mask, ((t == 1).astype(float) - (t == 2)) / s2, 0.0)]
        chi = np.where(~ok3, 0, np.where(r3 == 1, -1, 1)).astype(np.float64)
        out['E37m3'] = (lab, 6, np.array([0, 0, 0, 1, 1, 1]), mask, np.where(mask, t, 0), 3, F, chi)
    return out

def pair_sums(Fa, Fb):
    nF = len(Fa)
    return dict(n=len(Fa[0]) if nF else 0, sf=np.array([f.sum() for f in Fa]), sg=np.array([g.sum() for g in Fb]),
                sfg=np.array([[np.dot(f, g) for g in Fb] for f in Fa]), sff=np.array([np.dot(f, f) for f in Fa]),
                sgg=np.array([np.dot(g, g) for g in Fb]))

def do_segment(L, U, X1, cases):
    Ue = min(U + EXT, X1)
    seg = segment(L, Ue, cases=cases)
    p = seg['p']; nmain = int(np.searchsorted(p, U)); rec = {'L': L, 'U': U, 'nprimes': nmain}
    gaps = np.diff(p); gi = np.minimum(gaps // 2, GCAP)
    for name, (lab, A, alpha, mask, fib, C, F, chi) in schemes(seg).items():
        pre = name + '__'; lab = lab.astype(np.int64); fib = np.asarray(fib, np.int64)
        good = lab >= 0
        for k in (1, 2, 3):
            i = np.arange(0, min(nmain, len(p) - k)); j = i + k; m = good[i] & good[j]
            rec[pre + 'X%d' % k] = np.bincount(lab[i[m]] * A + lab[j[m]], minlength=A * A).reshape(A, A)
            mm = m & mask[i] & mask[j]; ii, jj = i[mm], j[mm]
            ps = pair_sums([f[ii] for f in F], [f[jj] for f in F])
            for key, v in ps.items(): rec[pre + 'H%d_%s' % (k, key)] = v
            # f of p_n by alpha of p_(n+k), and f of p_(n+k) by alpha of p_n
            na = int(alpha.max()) + 1
            m1 = m & mask[i]; aj = alpha[lab[j[m1]]]
            rec[pre + 'B%d_fwd_n' % k] = np.bincount(aj, minlength=na)
            rec[pre + 'B%d_fwd_s' % k] = np.array([np.bincount(aj, weights=f[i[m1]], minlength=na) for f in F])
            m2 = m & mask[j]; ai = alpha[lab[i[m2]]]
            rec[pre + 'B%d_bwd_n' % k] = np.bincount(ai, minlength=na)
            rec[pre + 'B%d_bwd_s' % k] = np.array([np.bincount(ai, weights=f[j[m2]], minlength=na) for f in F])
            if k == 1:
                g = gi[i[m]]
                rec[pre + 'Xg'] = np.bincount((g * A + lab[i[m]]) * A + lab[j[m]], minlength=(GCAP + 1) * A * A).reshape(GCAP + 1, A, A)
                g = gi[ii]
                rec[pre + 'G1_n'] = np.bincount(g, minlength=GCAP + 1)
                rec[pre + 'G1_sfg'] = np.array([np.bincount(g, weights=f[ii] * f[jj], minlength=GCAP + 1) for f in F])
                rec[pre + 'G1_sf'] = np.array([np.bincount(g, weights=f[ii], minlength=GCAP + 1) for f in F])
                rec[pre + 'G1_sg'] = np.array([np.bincount(g, weights=f[jj], minlength=GCAP + 1) for f in F])
                rec[pre + 'G1_sff'] = np.array([np.bincount(g, weights=f[ii] ** 2, minlength=GCAP + 1) for f in F])
                rec[pre + 'G1_sgg'] = np.array([np.bincount(g, weights=f[jj] ** 2, minlength=GCAP + 1) for f in F])
        i = np.arange(0, min(nmain, len(p) - 2)); m = good[i] & good[i + 1] & good[i + 2]
        rec[pre + 'X3t'] = np.bincount((lab[i[m]] * A + lab[i[m] + 1]) * A + lab[i[m] + 2], minlength=A ** 3).reshape(A, A, A)
        # mask subsequence
        mi = np.nonzero(mask)[0]; nm = int(np.searchsorted(mi, nmain))
        for k in (1, 2, 3):
            t = np.arange(0, min(nm, len(mi) - k)); ii, jj = mi[t], mi[t + k]
            ps = pair_sums([f[ii] for f in F], [f[jj] for f in F])
            for key, v in ps.items(): rec[pre + 'M%d_%s' % (k, key)] = v
        # runs of the fibre label inside [L, U) in the mask subsequence (first and last run dropped)
        s = fib[mi[:nm]]; hist = np.zeros((C, RUNCAP + 1), np.int64)
        if len(s) > 2:
            br = np.nonzero(np.diff(s) != 0)[0] + 1; starts = np.concatenate([[0], br]); ends = np.concatenate([br, [len(s)]])
            ln = (ends - starts)[1:-1]; val = s[starts][1:-1]
            np.add.at(hist, (val, np.minimum(ln, RUNCAP)), 1)
        rec[pre + 'runs'] = hist
        # block dispersion (primes in [L, U) only)
        inm = np.arange(nmain)
        for H in HS:
            nb = (U - L) // H; blk = (p[:nmain] - L) // H
            res = []
            for f in F + [chi]:
                T = np.bincount(blk, weights=f[:nmain], minlength=nb); V = np.bincount(blk, weights=f[:nmain] ** 2, minlength=nb)
                ok = V > 0; r = T[ok] ** 2 / V[ok]; T2 = T * T
                res.append([ok.sum(), r.sum(), (r * r).sum(), T2.sum(), V.sum(), (T2 * T2).sum(), (V * V).sum(), (T2 * V).sum(), nb])
            rec[pre + 'D%d' % H] = np.array(res, dtype=np.float64)
        # per-block (1e4) fibre transition tables of consecutive mask-mask pairs and abelian sign tables
        H = 10 ** 4; nb = (U - L) // H
        i = np.arange(0, min(nmain, len(p) - 1)); j = i + 1; m = mask[i] & mask[j]
        b = (p[i[m]] - L) // H
        rec[pre + 'W4'] = np.bincount((b * C + fib[i[m]]) * C + fib[j[m]], minlength=nb * C * C).reshape(nb, C, C).astype(np.int32)
        if name in ('C23', 'C31', 'C47', 'C39', 'G', 'E37m2'):
            m = good[i] & good[j] & (chi[i] != 0) & (chi[j] != 0); b = (p[i[m]] - L) // H
            si = (chi[i[m]] > 0).astype(np.int64); sj = (chi[j[m]] > 0).astype(np.int64)
            rec[pre + 'W4a'] = np.bincount((b * 2 + si) * 2 + sj, minlength=nb * 4).reshape(nb, 2, 2).astype(np.int32)
    return rec

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('x0', type=float); ap.add_argument('x1', type=float)
    ap.add_argument('--seg', type=float, default=5e6); ap.add_argument('--cases', default='C23,C31,C47,C39,G,E37')
    ap.add_argument('--max-seconds', type=float, default=540)
    a = ap.parse_args(); X0, X1, S = int(a.x0), int(a.x1), int(a.seg)
    assert X0 % S == 0 and X1 % S == 0 and S % 10 ** 6 == 0
    cases = tuple(a.cases.split(',')); os.makedirs(a.outdir, exist_ok=True); t0 = time.time()
    for L in range(X0, X1, S):
        fn = os.path.join(a.outdir, 'seg_%011d.npz' % L)
        if os.path.exists(fn): continue
        if time.time() - t0 > a.max_seconds: print('time budget reached; rerun to resume', flush=True); return
        ts = time.time(); rec = do_segment(L, L + S, X1, cases)
        np.savez_compressed(fn + '.tmp.npz', **{k: np.asarray(v) for k, v in rec.items()}); os.replace(fn + '.tmp.npz', fn)
        print('segment [%d, %d) %d primes  %.1fs' % (L, L + S, rec['nprimes'], time.time() - ts), flush=True)
    print('all segments done', flush=True)

if __name__ == '__main__':
    main()
