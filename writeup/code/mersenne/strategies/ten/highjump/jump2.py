"""The 85 record gaps (OEIS A002386 / A005250 b-files) against Wolf's formula (x / pi(x)) (2 ln pi(x) - ln x + ln C2),
pi(x) replaced by li(x) (relative error below 10^-4 for x > 10^6), Cramer's random-model maximum ln x ln pi(x), and
ln^2 x; all at x = p.  Also the largest g / ln^2 p.   usage: jump2.py out.txt"""
import sys, math
from scipy.special import expi
out = open(sys.argv[1], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
P = [int(l.split()[1]) for l in open('oeis/b002386.txt') if l.strip() and not l.startswith('#')]
G = [int(l.split()[1]) for l in open('oeis/b005250.txt') if l.strip() and not l.startswith('#')]
C2 = 0.6601618158468696
li = lambda x: float(expi(math.log(x)) - expi(math.log(2)))
r_w = []; r_c = []; r_l = []
for p, g in zip(P, G):
    if p < 10 ** 4: continue
    pi = li(p); w = p / pi * (2 * math.log(pi) - math.log(p) + math.log(C2)); c = math.log(p) * math.log(pi)
    r_w.append(g / w); r_c.append(g / c); r_l.append(g / math.log(p) ** 2)
n = len(r_w)
pr('%d records with p >= 10^4 (up to p = %d, gap %d)' % (n, P[-1], G[-1]))
pr('record / Wolf:              mean %.3f, min %.3f, max %.3f' % (sum(r_w) / n, min(r_w), max(r_w)))
pr('record / (ln p ln pi(p)):   mean %.3f, min %.3f, max %.3f' % (sum(r_c) / n, min(r_c), max(r_c)))
pr('record / ln^2 p:            mean %.3f, min %.3f, max %.3f' % (sum(r_l) / n, min(r_l), max(r_l)))
k = max(range(len(P)), key=lambda i: G[i] / math.log(P[i]) ** 2 if P[i] > 100 else 0)
pr('largest g / ln^2 p for p > 100: %.4f at p = %d (gap %d)' % (G[k] / math.log(P[k]) ** 2, P[k], G[k]))
