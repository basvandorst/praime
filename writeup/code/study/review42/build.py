import sys, time, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mypred
t = time.time()
lnh, lnPhiD, lnPhiO = mypred.build_tables()
np.savez(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/tables.npz', lnh=lnh, lnPhiD=lnPhiD, lnPhiO=lnPhiO)
print('built', time.time() - t)
# sanity: GUE closed form
from scipy.special import sici
P = mypred.Pred((lnh, lnPhiD, lnPhiO))
ls = np.array([0.5, 1, 2, 4, 8.0])
r = P.r
def gue(l):
    si, ci = sici(2 * np.pi * l)
    return (np.log(2 * np.pi * l) + np.euler_gamma + 1 - np.cos(2 * np.pi * l) - ci) / np.pi ** 2 + l * (1 - 2 / np.pi * si)
# GUE alone: set X=1 by L huge
print('GUE closed', gue(ls))
print('L=1e9 zeta', P.sigma2(ls, 1e9, [], mode='zeta'))
print('PhiD(0)', P.s_lnPhiD(0.0), 'lnPhiO(0)', P.lnPhiO[0], 'lnh small', lnh[:3])
