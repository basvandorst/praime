"""every predicted cell of predictions_blind.json with its reading (decision it enters) or 'no reading'. usage: b4_cells.py <pred.json> > cells.txt"""
import sys, json, math
P = json.load(open(sys.argv[1]))
print('# key | class | quantity | predicted | alternatives | noise sd | reading')
for k, e in P['L'].items():
    p = e['pred']; a = e['alt']; s = e['sd']
    print('%s | %s | sharp r=1 | %+.2f%+.2fi | H0 0; H-noroot %+.2f; H-conj %+.2f%+.2fi; H-c1 %+.2f%+.2fi | along %.2f, across %.2f | D1, D2, D3, D5%s%s' % (
        k, e['cls'], p[0], p[1], a['H_noroot'][0], a['H_conj'][0], a['H_conj'][1], a['H_c1'][0], a['H_c1'][1], s['along'], s['across'],
        ', D6' if (e['eps'][0] < -0.5 and abs(e['eps'][1]) < 1e-9) else '', ', D7' if abs(e['eps'][1]) > 1e-9 else ''))
    print('%s | %s | Hann r=1 | %+.2f%+.2fi | - | along %.2f | no reading (reported)' % (k, e['cls'], e['pred_hann'][0], e['pred_hann'][1], s['along_hann']))
    for nm, l in e['lines'].items():
        print('%s | %s | line r=%s | %+.2f%+.2fi | H0 0 | along %.2f | D4' % (k, e['cls'], nm, l['pred'][0], l['pred'][1], s['line_along']))
for k, e in P['products'].items():
    for rd, v in e['readings'].items():
        tag = dict(M1='D8', M2='D9', M3='D9', R_artin='D10', half_phase='D11', product_phase='D11')[rd]
        print('%s | %s | %s | %+.2f | %s | %.2f | %s' % (k, e['cls'], rd, v['pred'], 'H-prod0 0' if rd in ('M1', 'M2') else ('H0 0' if rd in ('R_artin', 'half_phase') else '-'), v['sd'], tag))
        if 'pred_hann' in v:
            print('%s | %s | %s Hann | %+.2f | - | %.2f | no reading (reported)' % (k, e['cls'], rd, v['pred_hann'], v.get('sd_hann', float('nan'))))
        if rd == 'R_artin':
            for nm, l in v['lines'].items():
                print('%s | %s | Artin line r=%s | %+.2f%+.2fi | - | %.2f | no reading (reported)' % (k, e['cls'], nm, l['pred'][0], l['pred'][1], l['sd']))
for pr in P['pairs']:
    print('%s+%s | pair | product phase | 0 | - | sqrt(n1+n2) | D12' % tuple(pr['keys']))
