"""Check the linear-response tail machinery on the class counts (S1 statistic), where the answer is known: slope -1/lam per
unit ln H (vq = 1 - (ln(H/q) + ...)/lam).  Response of the normalized class vector to a prime at x: z(c_x) = sqrt(C)(e_c - 1/C)."""
import numpy as np, math
lam = np.mean([math.log(6e9 + i*1e8 + 5e7) for i in range(24)])
def sp(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)
PR = sp(10**7)
def S_tail(hmax, y):
    p = PR[PR > y].astype(float); Ct = math.exp(np.sum(np.log1p(-1/(p-1)**2)) - 1/(1e7*math.log(1e7)))
    f = np.ones(hmax+1)
    for pp in PR[(PR > y) & (PR <= hmax)].tolist(): f[pp::pp] *= (pp-1)/(pp-2)
    return f*Ct
HS = [10**3, 10**4, 10**5, 10**6]; Hmax = max(HS); h = np.arange(1, Hmax+1)
for y, Q in ((11, 2310), (13, 30030)):
    St = S_tail(Hmax, y)
    ok = np.ones(Q, bool)
    for p in PR[PR <= y].tolist(): ok[::p] = False
    U = np.flatnonzero(ok)
    for q in (3, 5, 7):
        red = [a for a in range(q) if math.gcd(a, q) == 1]; C = len(red)
        A = np.zeros(Q)
        for c in range(C):
            full = np.zeros(Q); full[U] = math.sqrt(C) * ((U % q == red[c]).astype(float) - 1/C)
            F = np.fft.rfft(full); A += np.fft.irfft(F*np.conj(F), n=Q)
        A /= len(U)
        pref = (Q/len(U)) / ((C-1)*lam)
        base = A[h % Q] * (St[1:] - 1)
        W = {H: 2*pref*np.sum((1 - h[:H-1]/H)*base[:H-1]) for H in HS}
        print('y=%d q=%d  S1 tail: %s  per decade: %s   (theory -1/lam per unit lnH = %.4f per decade)' % (y, q, ' '.join('%+.4f' % W[H] for H in HS), ' '.join('%+.4f' % (W[HS[i+1]] - W[HS[i]]) for i in range(3)), -math.log(10)/lam))
