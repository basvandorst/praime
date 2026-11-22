import math, sys
sys.argv = ['x', '10000']
exec(open('c5_windows.py').read().split('l7, l10 =')[0])
l7, l10 = math.log(1e7), math.log(1e10)
for d in (0.1, 0.01):
    o = offwin(d, l7, l10); print(d, 'off %.4f  D+off %.4f' % (o, closed(d) + o), flush=True)
