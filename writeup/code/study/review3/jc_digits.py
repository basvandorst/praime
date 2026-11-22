import json, math
import mpmath as mp
mp.mp.dps = 30
E = json.load(open('/home/bas/scratch/study-review3/enum.json'))
U = {int(h): [mp.mpf(x) for x in v] for h, v in E['U'].items() if int(h) in (6, 30)}
def share(h, lam): return sum((-1) ** j * U[h][j] / lam ** (j + 1) for j in range(len(U[h])))
def D(u): return share(30, u) - share(6, u)
def F(Ue, U0):  # (N30-N6) e^{-Ue}, density per prime: dt/ln t, t = e^u
    return mp.quad(lambda u: mp.e ** (u - Ue) * D(u) / u, [U0, Ue - 40, Ue - 10, Ue])
for x0 in (mp.mpf(10) ** 3, mp.mpf(10) ** 15):
    r = mp.findroot(lambda Ue: F(Ue, mp.log(x0)), mp.log(mp.mpf('1.7427e35')))
    print('x0 = %s: crossing x = %s' % (mp.nstr(x0, 3), mp.nstr(mp.e ** r, 12)))
# variant: the pair density for gap h at t uses 1/(ln t ln(t+h)) ... approximated by lam -> ln t + h/(2t): negligible; also
# variant: N counted by the larger prime instead of the smaller: shift t by h
print('literature interval: 1.74274357320e35 .. 1.74274357327e35')
# sensitivity: perturb U[30] by a relative 1e-9 (the precision level of the tail factors G_k)
base = mp.findroot(lambda Ue: F(Ue, mp.log(10 ** 3)), mp.log(mp.mpf('1.7427e35')))
U[30] = [x * (1 + mp.mpf('1e-9')) for x in U[30]]
pert = mp.findroot(lambda Ue: F(Ue, mp.log(10 ** 3)), mp.log(mp.mpf('1.7427e35')))
print('relative shift of x per 1e-9 relative change of the gap-30 constants: %s' % mp.nstr(mp.e ** (pert - base) - 1, 4))
U[30] = [x / (1 + mp.mpf('1e-9')) for x in U[30]]
# convention: lambda = ln of the midpoint t + h/2 instead of ln t
def share2(h, t): return share(h, mp.log(t + h / 2))
def F2(Ue, U0): return mp.quad(lambda u: mp.e ** (u - Ue) * (share2(30, mp.e ** u) - share2(6, mp.e ** u)) / u, [U0, Ue - 40, Ue - 10, Ue])
r = mp.findroot(lambda Ue: F2(Ue, mp.log(10 ** 3)), base)
print('with lambda = ln(t + h/2): x = %s' % mp.nstr(mp.e ** r, 12))
