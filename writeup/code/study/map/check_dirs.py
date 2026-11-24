"""compare map.py report JSON (no LR) with the review's dirs.json / dirs_old.json, group by group"""
import json, sys, math
R44 = '../../../../writeup/code/study/review44/'
for rep, ref in (('rep_new_nolr.json', 'dirs.json'), ('rep_old_nolr.json', 'dirs_old.json')):
    A = json.load(open(rep))['cells']; B = json.load(open(R44 + ref))
    print('==', rep, 'vs', ref)
    for cell_ref, cell in (('5_1', '5_1'), ('7_1', '7_1'), ('3_1', '3_1'), ('10_2', '5_2')):
        G = {g['name']: g for g in A[cell]['groups']}
        mx = 0; n = 0
        for lab, v in B[cell_ref].items():
            g = G[lab]; assert g['ndir'] == v['n']
            for i, H in enumerate(('10000', '100000', '1000000')):
                mx = max(mx, abs(g['W'][i][0] - v['W'][H][0]), abs(g['W'][i][1] - v['W'][H][1])); n += 2
            mx = max(mx, abs(g['slope'][0] - v['slope'][0])); n += 1
            # the review's slope error is std/sqrt(n); the jackknife SE of a mean is identical
            mx = max(mx, abs(g['slope'][1] - v['slope'][1])); n += 1
        print('  %s (mine %s): %d numbers, max |difference| = %.2e' % (cell_ref, cell, n, mx))
