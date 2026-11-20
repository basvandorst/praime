"""Pomerance's heuristic for the number of Carmichael numbers (and base-2 pseudoprimes) below x:
x^(1 - (1 + o(1)) lnlnln x / lnln x).  With the o(1) left out, against the OEIS counts (A055553, A055550) and our own
counts up to 10^9; and the growth factor per decade of the counts.   usage: pomerance.py"""
import math
PSP = {3: 3, 4: 22, 5: 78, 6: 245, 7: 750, 8: 2057, 9: 5597, 10: 14884, 11: 38975, 12: 101629, 13: 264239, 14: 687007, 15: 1801533, 16: 4744920, 17: 12604009, 18: 33763684, 19: 91210364}
CAR = {3: 1, 4: 7, 5: 16, 6: 43, 7: 105, 8: 255, 9: 646, 10: 1547, 11: 3605, 12: 8241, 13: 19279, 14: 44706, 15: 105212, 16: 246683, 17: 585355, 18: 1401644, 19: 3381806, 20: 8220777, 21: 20138200, 22: 49679870, 23: 123381982, 24: 308279939}
print(' n   x^(1 - lnlnln x / lnln x)   Carmichael   overshoot   growth per decade: pseudoprimes  Carmichael')
for n in range(6, 25):
    x = 10.0 ** n; L = math.log(x); h = x ** (1 - math.log(math.log(L)) / math.log(L))
    g1 = '%.2f' % (PSP[n] / PSP[n - 1]) if n in PSP and n - 1 in PSP else '-'
    print('%2d %26.3e %12d %11.0f %28s %11.2f' % (n, h, CAR[n], h / CAR[n], g1, CAR[n] / CAR[n - 1]))
