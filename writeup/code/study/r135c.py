"""Round 135, after the forty-third review: pools the outputs of r135b.py (registered windows; six fresh windows at
2e11, 4e11, 5.5e11, 7e11, 8e11, 1e12) -> ../../data/study/r135c.json. usage: r135c.py"""
import os, json, math
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, '..', '..', 'data', 'study')
J = {n: json.load(open(os.path.join(D, 'r135b_%s.json' % n))) for n in ('reg', 'freshA', 'freshB')}
def pool(key, names):
    m = sum(J[n][key]['mean'] for n in names) / len(names); e = math.sqrt(sum(J[n][key]['err'] ** 2 for n in names)) / len(names)
    law = None if J[names[0]][key]['law'] is None else sum(J[n][key]['law'] for n in names) / len(names)
    return {'mean': m, 'err': e, 'law': law, 'z': None if law is None else (m - law) / e}
out = {'fresh': {}, 'all': {}}
for key in J['reg']:
    if key == 'windows': continue
    out['fresh'][key] = pool(key, ('freshA', 'freshB')); out['all'][key] = pool(key, ('reg', 'freshA', 'freshB'))
f = out['fresh']; a = out['all']
for BL in (10 ** 4, 10 ** 5, 10 ** 6):
    print('H = %.0e' % BL)
    print('  classes, 9 windows:  ' + '  '.join('q=%d %.4f+-%.4f law %.4f z %+.1f' % (q, a['S1_%d_%d' % (q, BL)]['mean'], a['S1_%d_%d' % (q, BL)]['err'], a['S1_%d_%d' % (q, BL)]['law'], a['S1_%d_%d' % (q, BL)]['z']) for q in (3, 4, 5, 7, 210)))
    print('  column, 9 windows:   ' + '  '.join('(%d,%d) %.4f+-%.4f formula %.4f z %+.1f' % (q, k, a['col_%d_%d_%d' % (q, k, BL)]['mean'], a['col_%d_%d_%d' % (q, k, BL)]['err'], a['col_%d_%d_%d' % (q, k, BL)]['law'], a['col_%d_%d_%d' % (q, k, BL)]['z']) for q, k in ((3, 1), (3, 2), (5, 1), (7, 1), (10, 2))))
    print('  within, 6 fresh:     ' + '  '.join('(%d,%d) %.4f+-%.4f' % (q, k, f['win_%d_%d_%d' % (q, k, BL)]['mean'], f['win_%d_%d_%d' % (q, k, BL)]['err']) for q, k in ((3, 1), (3, 2), (5, 1), (7, 1), (10, 2))))
    print('  within, 9 windows:   ' + '  '.join('(%d,%d) %.4f+-%.4f' % (q, k, a['win_%d_%d_%d' % (q, k, BL)]['mean'], a['win_%d_%d_%d' % (q, k, BL)]['err']) for q, k in ((3, 1), (3, 2), (5, 1), (7, 1), (10, 2))))
json.dump(out, open(os.path.join(D, 'r135c.json'), 'w'), indent=0)
