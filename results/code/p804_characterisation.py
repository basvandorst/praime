# Characterisation: a d-digit M ending in 9 has M + rev(M) palindromic iff
#   first digit 2, every inner symmetric pair (m_j, m_{d-1-j}) sums to 0 or 11, and for odd d the middle digit is ...?
def pal(x): s = str(x); return s == s[::-1]
for L in range(2, 9):
    real = {M for M in range(10**(L - 1) + 9, 10**L, 10) if pal(M + int(str(M)[::-1]))}
    mids = sorted({str(M)[L // 2] for M in real}) if L % 2 else None
    def ok(M):
        s = str(M)
        if s[0] != '2' or s[-1] != '9': return False
        for j in range(1, L // 2):
            if int(s[j]) + int(s[L - 1 - j]) not in (0, 11): return False
        if L % 2 and s[L // 2] not in mids: return False
        return True
    pred = {M for M in range(10**(L - 1) + 9, 10**L, 10) if ok(M)}
    print(f"{L} digits: {len(real)} real, characterisation {'EXACT' if pred == real else 'differs'}" + (f", middle digit in {mids}" if mids else ""))
