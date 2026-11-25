import math, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6); pari.set_real_precision(40)
P = 30000
s = np.ones(P + 1, bool); s[:2] = False
for i in range(2, int(P ** .5) + 1):
    if s[i]: s[i * i::i] = False
pr = np.nonzero(s)[0]
for (NL, k, tag) in ((29, 2, 'W29'), (1, 24, 'F24'), (23, 2, 'W23'), (31, 2, 'W31')):
    pari('mf = mfinit([%d,%d],%d)' % (NL, k, 1 if NL == 1 else 0)); pari('F = mfeigenbasis(mf)[1]'); pari('v = mfcoefs(F, %d)' % P)
    pari('pol = mffields(mf)[1]'); pari('vt = variable(pol)'); pari('rts = polroots(pol)')
    for i in (1, 2):
        pari('zz = real(rts[%d])' % i)
        pari('W = vector(#pr0, j, 0)') if False else None
        lam = np.array([float(pari('subst(lift(v[%d]), vt, zz)' % (int(p) + 1))) / p ** ((k - 1) / 2) for p in pr])
        X = 1.0 / pr
        lt = np.log(1 - (lam ** 2 - 1) * X)
        bad = (NL % pr == 0)
        lt[bad] = np.log((1 - lam[bad] ** 2 * X[bad]) / (1 - X[bad]))
        cum = np.cumsum(lt); m = pr >= P / 10; w = 1.0 / pr[m]
        c = math.exp(((cum[m] + 1 / (np.sqrt(pr[m]) * np.log(pr[m]))) * w).sum() / w.sum())
        sgn = 1.0
        if np.any(1 - (lam ** 2 - 1) * X <= 0):
            # negative factor at p = 2: product of signed factors
            fac = 1 - (lam ** 2 - 1) * X; fac[bad] = (1 - lam[bad] ** 2 * X[bad]) / (1 - X[bad])
            sgn = np.prod(np.sign(fac)); lt = np.log(np.abs(fac)); cum = np.cumsum(lt)
            c = sgn * math.exp(((cum[m] + 1 / (np.sqrt(pr[m]) * np.log(pr[m]))) * w).sum() / w.sum())
        print('%s level %d weight %d root %d: lambda_2 = %+.4f  c (direct, P=%d) = %+.5f' % (tag, NL, k, i, lam[0], P, c))
