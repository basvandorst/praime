"""12*H(m) for 0 <= m <= M (Hurwitz class numbers, H(0) set to 0), by counting reduced forms (a,b,c), 4ac - b^2 = m,
|b| <= a <= c, b >= 0 if |b| = a or a = c; weight 1/2 for (a,0,a), 1/3 for (a,a,a). Own code; checked against PARI in validate_trace.py."""
import numpy as np


def hurwitz12(M):
    H = np.zeros(M + 1, np.int32)
    a = 1
    while 3 * a * a <= M:
        for b in range(0, a + 1):                      # c = a
            m = 4 * a * a - b * b
            if m <= M:
                H[m] += 6 if b == 0 else (4 if b == a else 12)
        st = 4 * a
        for b in range(0, a + 1):                      # c > a
            s = 4 * a * (a + 1) - b * b
            if s > M:
                continue
            H[s::st] += 12 if (b == 0 or b == a) else 24
        a += 1
    return H
