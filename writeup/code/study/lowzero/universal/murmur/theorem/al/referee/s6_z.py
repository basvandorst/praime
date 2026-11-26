import sys; sys.path.insert(0,'/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/referee_al')
from s3_level import *
from s4_lemmaT import wstats
L = level(10103); N = L['N']; ps = L['ps']; dim = L['dim']; y = ps / N
for name, data, f in (('t_q', L['tq'], Dq), ('t_N', L['tN'], DN)):
    m, se = wstats(ps, N, data / dim)
    mod, _ = wstats(ps, N, np.array([f(v) for v in y]))
    z = (m - mod) / se
    print('%s N=%d: z per window ' % (name, N) + ' '.join('%5.1f' % v for v in z) + '  chi2(16) = %.1f' % float(np.sum(z ** 2)))
    # wrong-scale alternative for t_q: s < 2 sqrt(y) (variable p/N instead of p/Q)
    if name == 't_q':
        alt = lambda v: Dq(v / 2)
        ma, _ = wstats(ps, N, np.array([alt(v) for v in y]))
        print('   alternative scale x = p/N: chi2(16) = %.0f' % float(np.sum(((m - ma) / se) ** 2)))
        mp, _ = wstats(ps, N, np.array([Dq(v) + 48 * v for v in y]))   # Eisenstein sign flipped (+24y)
        print('   Eisenstein term +24y instead of -24y: chi2(16) = %.0f' % float(np.sum(((m - mp) / se) ** 2)))
