import numpy as np, json, sys
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
name = sys.argv[1]; excl = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 and sys.argv[2] else []
dat = {int(k): v for k, v in json.load(open(SCR + f'data_{name}.json')).items()}
Ds = [D for D in dat if D not in excl]
g = lambda k: np.array([dat[D][k] for D in Ds])
A = g('k3'); m3f = g('m3fix'); M3f = g('M3fix'); tr = g('k3trim'); dl = g('dfix'); M2 = g('M2fix')
n = len(Ds)
def ms(x): return x.mean(0), x.std(0, ddof=1)/np.sqrt(n)
for lab, x in (('A - m3fix (Palm edge + range)', A - m3f), ('m3fix - M3fix (mean subtraction)', m3f - M3f), ('A - M3fix (total bias)', A - M3f), ('trim3 - A (repo edge check)', tr - A)):
    m, s = ms(x)
    print(f'{name} n={n} {lab}:')
    print('   mean x1e-4', np.round(m*1e4, 2)); print('   se   x1e-4', np.round(s*1e4, 2))
# GLS-combined single offset of the total bias across 16 lengths
x = A - M3f; m = x.mean(0); C = np.cov(x.T)/n; o = np.ones(16)
b = o @ np.linalg.solve(C, m)/(o @ np.linalg.solve(C, o)); sb = 1/np.sqrt(o @ np.linalg.solve(C, o))
print(f'{name}: combined total bias (GLS one offset) = {b*1e4:+.2f} +- {sb*1e4:.2f} x1e-4')
x = A - m3f; m = x.mean(0); C = np.cov(x.T)/n
b = o @ np.linalg.solve(C, m)/(o @ np.linalg.solve(C, o)); sb = 1/np.sqrt(o @ np.linalg.solve(C, o))
print(f'{name}: combined Palm-edge part = {b*1e4:+.2f} +- {sb*1e4:.2f} x1e-4')
print('sd of delta (fixed band) per set at l=1,4,8:', np.round(dl.std(0)[[1, 7, 15]], 5))
