"""L(1, Sym^2 f) for the level-1 eigenforms, built in ARITHMETIC normalisation with exact coefficients (GP direuler),
checked by lfuncheckfeq; then c = B(1)/L(1,Sym^2 f) with B over p <= 3000 from the constants files' lambda's (recomputed here).
Also the Petersson relation: L(1, Sym^2 f) against <f,f> (PARI mfpetersson)."""
import sys, os, math, cmath, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from ulib import *
pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6); pari.default('realprecision', 57)
P0 = 3000
for key, spec, T2, fam in CATALOGUE:
    if fam != 'GL2 level 1': continue
    if len(sys.argv) > 1 and key not in sys.argv[1:]: continue
    t0 = time.time()
    _, N, k, chi, idx, emb = spec
    pari('mfU = mfinit([1,%d],0); FU = mfeigenbasis(mfU)[%d]; AU = mfcoefs(FU, 60000)' % (k, idx + 1))
    if emb is None:
        pari('aU = AU')
    else:
        pari('aU = mfembed(FU, AU)[%d]' % (emb + 1)) if False else pari('EM = mfembed(FU, AU); aU = if(type(EM)=="t_MAT", EM[%d,], if(type(EM[1])=="t_VEC", EM[%d], EM))' % (emb + 1, emb + 1))
    pari('kU = %d' % k)
    pari('S2U = lfuncreate([n -> direuler(p=2, n, 1/((1-p^(kU-1)*X)*(1-(real(aU[p+1])^2-2*p^(kU-1))*X+p^(2*kU-2)*X^2))), 0, [2-kU,0,1], 2*kU-1, 1, 1])')
    feq = float(pari('lfuncheckfeq(S2U)'))
    LS1 = to_c(pari('lfun(S2U, kU)'))
    # Petersson norm (PARI convention), for the relation
    try:
        pet = to_c(pari('mfpetersson(mfsymbol(mfU, FU))')) if emb is None else None
    except Exception as e:
        pet = None
    # B over p <= P0 from exact a_p
    lnB = 0j
    lam = {}
    for p in primes_upto(P0):
        p = int(p)
        ap = float(pari('real(aU[%d])' % (p + 1)))
        l = ap / p ** ((k - 1) / 2); lam[p] = l
        X = 1 / p
        F = (1 - X) * (1 - (l * l - 1) * X)
        A = 1 / ((1 - X) * (1 - (l * l - 2) * X + X * X))
        lnB += cmath.log(F * A) - math.log(1 - X)
    B = cmath.exp(lnB).real
    c = B / LS1.real
    d = json.load(open('constants_%s.json' % key))
    rel = None
    if pet is not None:
        rel = LS1.real / pet.real * math.gamma(k) / (4 * math.pi) ** k
    print('%-5s k=%d feq=%.0f L(1,Sym2)=%.12f B=%.8f c=%.6f | direct smooth %.6f end %.6f | lambda_2 %.4f (file %.4f) | <f,f>=%s  L(1,Sym2)*Gamma(k)/((4pi)^k <f,f>) = %s  (%.0fs)' % (
        key, k, feq, LS1.real, B, c, d['c_direct_smooth'][0], d['c_direct_end'][0], lam[2], -d['poly_2_3']['2'][1][0],
        None if pet is None else '%.6e' % pet.real, None if rel is None else '%.10f' % rel, time.time() - t0), flush=True)
    d.update(c_exact=[c, 0.0], A1=LS1.real, B=[B, 0.0], A_how='GP direuler Sym2, arithmetic normalisation, exact a_p', A_feq=feq,
             c=c, c_method='exact: B(1)/L(1,Sym2 f) (arithmetic GP build), B over p <= %d' % P0, petersson=None if pet is None else pet.real,
             petersson_relation=rel)
    # rescale the lines to the new c
    eps = complex(*d['eps'])
    old = d.pop('lines_per_unit')
    d['lines_per_unit'] = old  # recomputed in s2_freeze from local data
    json.dump(d, open('constants_%s.json' % key, 'w'), indent=1)
