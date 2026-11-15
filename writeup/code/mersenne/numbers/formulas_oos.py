"""Follow-up on floor(c*r^n): the 'max OOS over tied in-sample fits' statistic, conditioned on the in-sample
best count (a smaller in-sample best leaves a larger feasible region for the far strips), and the honest OOS of a
single fit chosen without the later terms: the centre of the in-sample feasible region (mean of the tied vertices)."""
import numpy as np, math
src = open('formulas.py').read()
exec(src[src.index('def strips'):src.index('def exact_check')])
D = np.load('seqs.npz'); real = D['real']
def stats(a):
    lo, hi, phi = strips(a, 'F1')
    m = len(a); I, J = np.triu_indices(m, 1); U = []; V = []
    for bi in (lo, hi):
        for bj in (lo, hi):
            v = (bj[J] - bi[I]) / (phi[J] - phi[I]); U.append(bi[I] - phi[I] * v); V.append(v)
    U = np.concatenate(U); V = np.concatenate(V)
    best = (-1, 0, 0)
    pts = []
    for du in (-1, 1):
        for dv in (-1, 1):
            uu = U + du * 1e-13 * np.maximum(1, np.abs(U)); vv = V + dv * 1e-13 * np.maximum(1, np.abs(V)) / 60
            val = uu[:, None] + phi[None, :] * vv[:, None]
            ins = (val >= lo[None, :]) & (val < hi[None, :])
            pts.append((uu, vv, ins[:, :25].sum(1), ins[:, 25:].sum(1)))
    uu = np.concatenate([p[0] for p in pts]); vv = np.concatenate([p[1] for p in pts])
    ci = np.concatenate([p[2] for p in pts]); co = np.concatenate([p[3] for p in pts])
    mi = ci.max(); sel = ci == mi
    # honest fit: centre of the tied vertices of one optimal region (the most common in-sample hit set is not tracked;
    # use the vertex mean restricted to points with the in-sample optimum and check it is still optimal)
    uc, vc = uu[sel].mean(), vv[sel].mean()
    val = uc + phi * vc
    ins = (val >= lo) & (val < hi)
    return int(mi), int(co[sel].max()), int(ins[:25].sum()), int(ins[25:].sum())
r = stats(real)
print('real: in-sample best', r[0], 'max OOS over ties', r[1], '| centre fit: in-sample', r[2], 'OOS', r[3])
for nm in ('nullA', 'nullB'):
    S = np.array([stats(s) for s in D[nm]])
    cond = S[S[:, 0] <= r[0]]
    print(f'{nm}: P(max OOS over ties >= {r[1]}) = {(S[:, 1] >= r[1]).mean():.3f}; among nulls with in-sample best <= {r[0]} '
          f'({len(cond)} of {len(S)}): {(cond[:, 1] >= r[1]).mean():.3f}; centre-fit OOS mean {S[:, 3].mean():.3f}, max {S[:, 3].max()}')
