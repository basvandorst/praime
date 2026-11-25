"""consec_analyze.py: readings from the segment files of consec_tables.py.

usage: consec_analyze.py DIR OUT.json [--x0 X0 --x1 X1]   (only segments with X0 <= L < X1 are read)
Readings per scheme (all definitions in consec.md section 3):
  marg       label frequencies against Chebotarev / Sato-Tate
  ab_b       abelian correlation b_k = E[chi(p_n) chi(p_(n+k))] (labels with a real abelian character)
  corr       fibre-fibre interaction correlations r(f_i(p_n), f_j(p_(n+k))), k = 1, 2, 3 (consecutive in all primes), z = r sqrt(n)
  corrM      the same for the k-th next fibre prime
  indep      contingency chi^2 of the fibre labels inside each abelian fibre pair (s, t), k = 1, 2, 3
  nbr        f(p_n) against the abelian datum of p_(n+k) (and f(p_(n+k)) against that of p_n): contrast z
  gap        lag-1 interaction by gap g = 2, 4, ...: chi^2 = sum z_g^2 over gaps with n_g >= 400; and the contrast
             gaps = 0 mod 2q (q the abelian modulus) against the rest
  triple     first-third interaction given the middle label (discrete schemes), chi^2
  runs       run lengths of the fibre label in the fibre subsequence against the geometric law
  w4         per-block (1e4) within / column parts of the fibre transition table (r135b method) and of the abelian sign table
  disp       block dispersion R_H = sum T^2 / sum V per function and for the abelian sign, H = 1e2 .. 1e6, per height bin
"""
import sys, os, json, math, glob, argparse
import numpy as np
from scipy.stats import chi2 as CHI2

CHEB_COND = {   # Chebotarev / Sato-Tate law of the fibre class code inside the fibre
    'C23': [1 / 3, 2 / 3], 'C31': [1 / 3, 2 / 3], 'C47': [1 / 5, 2 / 5, 2 / 5], 'C39': [1 / 2, 1 / 2],
    'G': [1 / 4] * 4, 'E37s': [1 / 2, 1 / 2], 'E37m2': [1 / 3, 2 / 3]}
ALPHA = {'C23': [0, 1, 1], 'C31': [0, 1, 1], 'C47': [0, 1, 1, 1], 'C39': [0, 1, 1, 2], 'G': [0, 1, 1, 1, 1],
         'E37s': [0, 0, 0], 'E37m2': [0, 1, 1], 'E37m3': [0, 0, 0, 1, 1, 1]}
FIBRE = {'C23': [1, 2], 'C31': [1, 2], 'C47': [1, 2, 3], 'C39': [1, 2], 'G': [1, 2, 3, 4], 'E37s': [0, 1], 'E37m2': [1, 2],
         'E37m3': None}
# abelian sign chi per label (+1 / -1), modulus q of the abelian character
CHI = {'C23': ([-1, 1, 1], 23), 'C31': ([-1, 1, 1], 31), 'C47': ([-1, 1, 1, 1], 47), 'G': ([-1, 1, 1, 1, 1], 4),
       'E37m2': ([-1, 1, 1], 37), 'E37m3': ([1, 1, 1, -1, -1, -1], 3), 'C39': ([-1, 1, 1, 1], 39)}
MARG = {'C23': [1 / 2, 1 / 6, 1 / 3], 'C31': [1 / 2, 1 / 6, 1 / 3], 'C47': [1 / 2, 1 / 10, 1 / 5, 1 / 5],
        'C39': [1 / 2, 1 / 8, 1 / 8, 1 / 4], 'G': [1 / 2] + [1 / 8] * 4, 'E37m2': [1 / 2, 1 / 6, 1 / 3],
        'E37m3': [1 / 8, 3 / 16, 3 / 16, 1 / 4, 1 / 8, 1 / 8]}
GCAP = 100; HS = (100, 1000, 10000, 100000, 1000000)

def load(d, x0, x1):
    fs = sorted(glob.glob(os.path.join(d, 'seg_*.npz'))); segs = []
    for f in fs:
        z = np.load(f); L = int(z['L'])
        if x0 <= L < x1: segs.append({k: z[k] for k in z.files})
    return segs

def ssum(segs, key): return sum(s[key].astype(np.float64) for s in segs)

def corr_from(n, sf, sg, sfg, sff, sgg, centred=True):
    if n < 2: return float('nan'), float('nan')
    if centred:
        c = sfg - sf * sg / n; vf = sff - sf * sf / n; vg = sgg - sg * sg / n
    else: c, vf, vg = sfg, sff, sgg
    r = c / math.sqrt(vf * vg); return float(r), float(r * math.sqrt(n))

def contingency(T):
    T = np.asarray(T, float); T = T[T.sum(1) > 0][:, T.sum(0) > 0]
    if T.shape[0] < 2 or T.shape[1] < 2: return 0.0, 0
    E = T.sum(1, keepdims=True) * T.sum(0, keepdims=True) / T.sum()
    return float(((T - E) ** 2 / E).sum()), (T.shape[0] - 1) * (T.shape[1] - 1)

def within_column(W, C):
    """r135b method on per-block tables W (NB, C, C): pooled P, chi^2 per block, column part, within part."""
    W = W.astype(np.float64); NB = W.shape[0]
    Na = W.sum(2); Ph = W.sum(0) / W.sum(0).sum(1, keepdims=True); E = Na[:, :, None] * Ph[None]
    ok = (E > 0).all((1, 2)); W, E, Na = W[ok], E[ok], Na[ok]; NB = W.shape[0]
    chi = ((W - E) ** 2 / E).sum((1, 2)); cv = (W - E).sum(1)
    M = np.einsum('ba,aij->bij', Na, np.array([np.diag(Ph[i]) - np.outer(Ph[i], Ph[i]) for i in range(C)]))
    col = np.einsum('bi,bij,bj->b', cv, np.linalg.pinv(M), cv); f = NB / (NB - 1.0)
    colv = col / (C - 1) * f; win = (chi - col) / (C - 1) ** 2 * f
    return dict(blocks=int(NB), col=float(colv.mean()), col_err=float(colv.std() / math.sqrt(NB)),
                within=float(win.mean()), within_err=float(win.std() / math.sqrt(NB)), P=Ph.tolist())

def ratio_err(D):
    """R = sum T^2 / sum V with delta-method error from block moments (blocks independent)."""
    nb = D[8]; ST2, SV, ST4, SV2, ST2V = D[3], D[4], D[5], D[6], D[7]
    if SV <= 0: return float('nan'), float('nan')
    R = ST2 / SV; mT, mV = ST2 / nb, SV / nb
    vT = ST4 / nb - mT ** 2; vV = SV2 / nb - mV ** 2; cTV = ST2V / nb - mT * mV
    var = (vT - 2 * R * cTV + R * R * vV) / (nb * mV * mV)
    return float(R), float(math.sqrt(max(var, 0)))

def analyse(segs):
    out = {'segments': len(segs), 'range': [int(min(s['L'] for s in segs)), int(max(s['U'] for s in segs))],
           'nprimes': int(sum(int(s['nprimes']) for s in segs))}
    names = sorted({k.split('__')[0] for k in segs[0] if '__' in k})
    for nm in names:
        pre = nm + '__'; R = {}
        X = {k: ssum(segs, pre + 'X%d' % k) for k in (1, 2, 3)}; A = X[1].shape[0]
        # marginals (first index of X1)
        cnt = X[1].sum(1); N = cnt.sum()
        if nm in MARG:
            P = np.array(MARG[nm]); R['marg'] = dict(counts=cnt.tolist(), cheb=P.tolist(),
                                                    z=((cnt - N * P) / np.sqrt(N * P * (1 - P))).tolist())
        else:
            R['marg'] = dict(counts=cnt.tolist())
        # abelian correlation
        if nm in CHI:
            ch = np.array(CHI[nm][0], float)
            R['ab_b'] = {str(k): dict(b=float(ch @ X[k] @ ch / X[k].sum()), n=float(X[k].sum()),
                                      mean_chi=float((ch @ X[k].sum(1)) / X[k].sum())) for k in (1, 2, 3)}
        # interaction correlations (centred and uncentred)
        nF = segs[0][pre + 'H1_sf'].shape[0]
        for tag in ('H', 'M'):
            res = {}
            for k in (1, 2, 3):
                n = ssum(segs, pre + '%s%d_n' % (tag, k)); sf = ssum(segs, pre + '%s%d_sf' % (tag, k)); sg = ssum(segs, pre + '%s%d_sg' % (tag, k))
                sfg = ssum(segs, pre + '%s%d_sfg' % (tag, k)); sff = ssum(segs, pre + '%s%d_sff' % (tag, k)); sgg = ssum(segs, pre + '%s%d_sgg' % (tag, k))
                for i in range(nF):
                    for j in range(nF):
                        r, z = corr_from(float(n), sf[i], sg[j], sfg[i, j], sff[i], sgg[j])
                        r0, z0 = corr_from(float(n), sf[i], sg[j], sfg[i, j], sff[i], sgg[j], centred=False)
                        res['k%d_f%d_f%d' % (k, i, j)] = dict(n=float(n), r=r, z=z, r_uncentred=r0, z_uncentred=z0,
                                                             mean_f=float(sf[i] / n), mean_g=float(sg[j] / n))
            R['corr' if tag == 'H' else 'corrM'] = res
        # contingency inside abelian fibre pairs (discrete schemes)
        if FIBRE.get(nm) is not None or nm == 'E37m3':
            al = np.array(ALPHA[nm]); ind = {}
            for k in (1, 2, 3):
                tot, dof = 0.0, 0; parts = {}
                for s in sorted(set(al.tolist())):
                    for t in sorted(set(al.tolist())):
                        rows = np.nonzero(al == s)[0]; cols = np.nonzero(al == t)[0]
                        if nm == 'E37s': rows = cols = np.array([0, 1])
                        c2, d = contingency(X[k][np.ix_(rows, cols)])
                        if d: parts['%d,%d' % (s, t)] = [c2, d]; tot += c2; dof += d
                ind[str(k)] = dict(chi2=tot, dof=dof, p=float(CHI2.sf(tot, dof)) if dof else None, parts=parts)
            R['indep'] = ind
        # neighbour abelian datum
        nbr = {}
        for k in (1, 2, 3):
            for side in ('fwd', 'bwd'):
                nn = ssum(segs, pre + 'B%d_%s_n' % (k, side)); ss = ssum(segs, pre + 'B%d_%s_s' % (k, side))
                if len(nn) < 2: continue
                for i in range(nF):
                    m = ss[i] / np.maximum(nn, 1)
                    # contrast of the first abelian value against the rest (f has variance ~1 inside the fibre)
                    for a in range(1, len(nn)):
                        z = (m[a] - m[0]) / math.sqrt(1 / max(nn[a], 1) + 1 / max(nn[0], 1))
                        nbr['k%d_%s_f%d_a%d' % (k, side, i, a)] = dict(means=m.tolist(), n=nn.tolist(), z=float(z))
        R['nbr'] = nbr
        # gap-conditioned interaction (function 0 and all functions on the diagonal)
        gn = ssum(segs, pre + 'G1_n'); gfg = ssum(segs, pre + 'G1_sfg'); gf = ssum(segs, pre + 'G1_sf'); gg = ssum(segs, pre + 'G1_sg')
        gff = ssum(segs, pre + 'G1_sff'); ggg = ssum(segs, pre + 'G1_sgg'); gap = {}
        q = CHI[nm][1] if nm in CHI else None
        for i in range(nF):
            zs, used = [], []
            for g in range(1, GCAP):
                if gn[g] < 400: continue
                r, z = corr_from(float(gn[g]), gf[i, g], gg[i, g], gfg[i, g], gff[i, g], ggg[i, g]); zs.append(z); used.append(2 * g)
            zs = np.array(zs); ent = dict(gaps=used, z=zs.tolist(), chi2=float((zs ** 2).sum()), dof=len(zs),
                                         p=float(CHI2.sf((zs ** 2).sum(), len(zs))), max_abs_z=float(np.abs(zs).max()))
            if q:
                qq = 2 * q // math.gcd(2, q) if q % 2 else q   # gaps are even: g = 0 mod q  <=>  g = 0 mod lcm(2, q)
                sel = np.array([(2 * g) % qq == 0 for g in range(GCAP + 1)]); sel[0] = False; sel[GCAP] = False
                for lab_, msk in (('mult', sel), ('rest', ~sel & (np.arange(GCAP + 1) > 0) & (np.arange(GCAP + 1) < GCAP))):
                    n = gn[msk].sum(); r, z = corr_from(float(n), gf[i, msk].sum(), gg[i, msk].sum(), gfg[i, msk].sum(), gff[i, msk].sum(), ggg[i, msk].sum())
                    ent['gap_' + lab_] = dict(n=float(n), r=r, z=z)
            gap['f%d' % i] = ent
        R['gap'] = gap
        # triples (discrete): interaction of first and third fibre labels given the middle label
        if FIBRE.get(nm) is not None:
            T3 = ssum(segs, pre + 'X3t'); fb = FIBRE[nm]; tr = {}
            for mlab in range(A):
                c2, d = contingency(T3[np.ix_(fb, [mlab], fb)][:, 0, :])
                tr[str(mlab)] = [c2, d]
            tot = sum(v[0] for v in tr.values()); dof = sum(v[1] for v in tr.values())
            R['triple'] = dict(parts=tr, chi2=tot, dof=dof, p=float(CHI2.sf(tot, dof)))
        # runs
        if nm in CHEB_COND:
            H = ssum(segs, pre + 'runs'); C = H.shape[0]; res = {}; tot, dof = 0.0, 0
            for c in range(C):
                h = H[c]; nr = h.sum()
                if nr == 0: continue
                pc = CHEB_COND[nm][c]
                ell = np.arange(1, 41); pr = pc ** (ell - 1) * (1 - pc); pr[-1] = pc ** 39   # last bin: >= 40
                obs = h[1:41]; exp = nr * pr
                # pool the tail so every expected count >= 20
                cut = int(np.argmax(exp < 20)) if (exp < 20).any() else len(exp)
                o = np.concatenate([obs[:cut], [obs[cut:].sum()]]); e = np.concatenate([exp[:cut], [exp[cut:].sum()]])
                c2 = float(((o - e) ** 2 / e).sum()); d = len(o) - 1
                res[str(c)] = dict(chi2=c2, dof=d, mean_len=float((h[1:] * np.arange(1, 41)).sum() / nr), geo_mean=1 / (1 - pc))
                tot += c2; dof += d
            res['total'] = dict(chi2=tot, dof=dof, p=float(CHI2.sf(tot, dof)))
            R['runs'] = res
        # per-block within / column
        W = np.concatenate([s[pre + 'W4'] for s in segs]); C = W.shape[1]
        R['w4'] = within_column(W, C)
        if pre + 'W4a' in segs[0]:
            Wa = np.concatenate([s[pre + 'W4a'] for s in segs]); R['w4a'] = within_column(Wa, 2)
        # dispersion per height bin (decades of segments) and overall
        Ls = np.array([int(s['L']) for s in segs]); disp = {}
        edges = sorted(set([out['range'][0]] + [e for e in (10 ** 7, 3 * 10 ** 7, 10 ** 8, 3 * 10 ** 8, 10 ** 9) if out['range'][0] < e < out['range'][1]] + [out['range'][1]]))
        for H in HS:
            for b0, b1 in [(out['range'][0], out['range'][1])] + list(zip(edges[:-1], edges[1:])):
                sel = [s for s, L in zip(segs, Ls) if b0 <= L < b1]
                if not sel: continue
                D = ssum(sel, pre + 'D%d' % H); ent = {}
                for i in range(D.shape[0]):
                    Rv, Re = ratio_err(D[i]); nb = D[i, 0]
                    ent['f%d' % i if i < D.shape[0] - 1 else 'chi'] = dict(R=Rv, err=Re, mean_r=float(D[i, 1] / max(nb, 1)),
                                                                        mean_r_err=float(math.sqrt(max(D[i, 2] / max(nb, 1) - (D[i, 1] / max(nb, 1)) ** 2, 0) / max(nb, 1))), blocks=float(nb))
                disp['H%d_%d_%d' % (H, b0, b1)] = ent
        R['disp'] = disp
        out[nm] = R
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('dir'); ap.add_argument('out'); ap.add_argument('--x0', type=float, default=0); ap.add_argument('--x1', type=float, default=1e30)
    a = ap.parse_args(); segs = load(a.dir, int(a.x0), int(a.x1))
    res = analyse(segs); json.dump(res, open(a.out, 'w'), indent=1)
    print('segments', res['segments'], 'range', res['range'], 'primes', res['nprimes'])
