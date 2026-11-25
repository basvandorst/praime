import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import math, numpy as np
from g2data import load
rows = load()
w = {r['key']: r['c'] ** 2 for r in rows}; W = sum(w.values())
pos = sorted([r for r in rows if r['eps'] > 0], key=lambda r: -r['c'])
neg = sorted([r for r in rows if r['eps'] < 0], key=lambda r: -r['c'])
print('weight share eps+1: %.3f' % (sum(w[r['key']] for r in pos) / W))
for r in pos: print('+ %-14s %4s c %.3f w %.3f G/pred %.3f' % (r['key'], r['rnd'], r['c'], w[r['key']] / W, r['G'].real / (-r['eps'] * r['c'] * 150 / 2 / math.pi)))
# eps -1: every 6th by c (spread), plus the two most negative c
sel = neg[::6] + neg[-2:]
print('eps-1 picks:', ' '.join(r['key'] for r in sel))
print('eps-1 picks c:', ' '.join('%.2f' % r['c'] for r in sel))
open('sample.txt', 'w').write(' '.join(r['key'] for r in pos + sel))
