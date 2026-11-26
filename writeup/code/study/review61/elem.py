# review 162, task 3b: elementary classes of primes, per prime (own root counts on a random subset of fields).
# H1: 4p = t^2 - D n^2 (principal form of the order of conductor n of the quadratic resolvent Q(sqrt D)); thresholds y = n^2/4
#     (D < 0: p >= |D| n^2/4; D > 0: p <= D n^2/4). Such p split completely in the ring class field, so (if D_R = D_K) in the Galois
#     closure of the cubic resolvent: Frobenius in V4, a in {3, -1}. In the S3 cubic case the same primes have a = 2 always.
# H3: p = |g(t)| for integers t (g = polredabs polynomial): a principal degree-1 prime, a >= 0.
# Also: distribution of a in the windows just right of u = 1/3, 2/3 (cube-line spikes) against neighbouring windows.
import sys, math, time, numpy as np
sys.argv = [sys.argv[0], '1']
exec(open('/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/roots.py').read().split("rows = [l.split()")[0])
NF = {1: 800, 2: 700}; t0 = time.time(); res = {}
for r2 in (1, 2):
    rows = [l.split() for l in open(R + 'fields_r%d.txt' % r2)]; sel = np.sort(np.random.default_rng(99 + r2).choice(len(rows), NF[r2], replace=False))
    acc = {'all': [0, 0.0], 'H3': [0, 0.0]}; hist = {}; Hn = {n: [0, 0.0, 0.0] for n in range(1, 6)}   # count, sum a, sum a over the same y-window for all primes
    fine_all = np.zeros((2, 384)); fine_noH1 = np.zeros((2, 384)); fine_H1 = np.zeros((2, 384)); win = {}
    for i in sel:
        D = int(rows[i][0]); g = pari(rows[i][1]); A = abs(D); cf = coeffs(g)
        p = PR[PR <= 8 * A]; p = p[(A % p != 0) & (p >= 5)]
        ind = math.isqrt(abs(int(pari.poldisc(g))) // A); p = p[ind % p != 0]          # drop the few index primes (not needed here)
        a = nroots_trace(cf, p) - 1
        h1 = np.zeros(len(p), bool)
        for n in range(1, 6):
            s = 4 * p - A * n * n if D < 0 else A * n * n - 4 * p
            ok = s >= 0; r = np.zeros(len(p), np.int64); r[ok] = np.round(np.sqrt(s[ok].astype(float))).astype(np.int64)
            sq = ok & (r * r == s); h1 |= sq
            Hn[n][0] += int(sq.sum()); Hn[n][1] += float(a[sq].sum())
        T = int((8 * A)**0.25) + 3; vals = set(abs(int(pari.subst(g, 'x', t))) for t in range(-T, T + 1)); h3 = np.isin(p, np.array(sorted(vals), np.int64))
        acc['all'][0] += len(p); acc['all'][1] += float(a.sum()); acc['H3'][0] += int(h3.sum()); acc['H3'][1] += float(a[h3].sum())
        y = p / A; u = np.cbrt(y); j = np.minimum((u * 192).astype(int), 383); wl = np.log(p.astype(float)); sA = math.sqrt(A)
        for arr, msk in ((fine_all, np.ones(len(p), bool)), (fine_noH1, ~h1), (fine_H1, h1)):
            arr[0] += np.bincount(j[msk], weights=(wl * sA * a)[msk], minlength=384); arr[1] += np.bincount(j[msk], weights=wl[msk], minlength=384)
        for lab, (u0, u1) in (('after 1/3', (1 / 3, 1 / 3 + 1 / 192)), ('control 1/3', (1 / 3 + 2 / 192, 1 / 3 + 6 / 192)), ('after 2/3', (2 / 3, 2 / 3 + 1 / 192)),
                              ('control 2/3', (2 / 3 + 2 / 192, 2 / 3 + 6 / 192))):
            mk = (u > u0) & (u <= u1); hh = win.setdefault(lab, np.zeros(5)); hh += np.bincount(a[mk] + 1, minlength=5)[:5]
    print('r2=%d: %d fields, %d primes, mean a %.4f; H3 (values of g) %d primes, mean a %.3f; %.0f s' % (r2, len(sel), acc['all'][0], acc['all'][1] / acc['all'][0],
          acc['H3'][0], acc['H3'][1] / max(1, acc['H3'][0]), time.time() - t0), flush=True)
    print('   H1 (4p = t^2 - D n^2): ' + '; '.join('n=%d: %d primes, mean a %.3f' % (n, v[0], v[1] / max(1, v[0])) for n, v in Hn.items()))
    for lab, hh in win.items():
        tot = hh.sum(); print('   %-12s %6d primes: fraction a=-1 %.3f, a=0 %.3f, a=1 %.3f, a=3 %.4f; mean a %.4f' % (lab, tot, hh[0] / tot, hh[1] / tot, hh[2] / tot, hh[4] / tot, (hh * np.arange(-1, 4)).sum() / tot))
    m_all = fine_all[0] / fine_all[1]; m_no = fine_noH1[0] / fine_noH1[1]
    for y0 in (0.25, 1.0, 2.25, 4.0, 6.25):
        k0 = int(y0**(1 / 3) * 192); ks = range(k0 - 2, k0 + 3)
        print('   y=%.2f fine bins u=%s: all %s | without H1 %s' % (y0, ' '.join('%.3f' % ((k + 0.5) / 192) for k in ks), ' '.join('%6.2f' % m_all[k] for k in ks), ' '.join('%6.2f' % m_no[k] for k in ks)))
    np.savez(R + 'elem_r%d.npz' % r2, fine_all=fine_all, fine_noH1=fine_noH1, fine_H1=fine_H1)
