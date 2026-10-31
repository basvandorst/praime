# Puzzle 93: if N's base-b digits, reversed, are its base-c digits (k digits, leading zero allowed on the
# reversed side), then N * (b^(k-1) - 1) = 0 mod (bc - 1), so N = 0 mod (bc-1)/gcd(bc-1, b^(k-1)-1).
from math import gcd

def digits(n, base):
    d = []
    while n:
        d.append(n % base); n //= base
    return d[::-1] or [0]

def is_reversal(N, b, c):
    db = digits(N, b)
    s = db[::-1]                     # reversed base-b digits, read in base c
    if max(s) >= c: return False
    return sum(x * c**(len(s) - 1 - i) for i, x in enumerate(s)) == N

def forced(N, b, c):
    k = len(digits(N, b))
    m = b * c - 1
    return m // gcd(m, b**(k - 1) - 1)

# 1. the page's examples (decimal b = 10, other base c); 774 is printed as 2*3^2*32 on the page, it is 2*3^2*43
page = {16: [53, 371, 5141, 99481, 8520280], 4: [13], 7: [23, 46, 2116, 15226], 8: [1527465],
        9: [445, 313725], 12: [315231], 13: [43, 86, 774], 14: [834], 19: [21, 42, 63, 84, 441, 882, 1540, 3290, 7721, 9471]}
for c, Ns in page.items():
    for N in Ns:
        assert is_reversal(N, 10, c), (N, c)
        f = forced(N, 10, c)
        assert N % f == 0, (N, c, f)
        print(f"base {c:2d}: {N:>8d} reversal, forced divisor {f:>4d}, N/{f} = {N // f}")

# 2. all decimal/hex reversals with more than one digit: decimal digits <= 9 means hex digits <= 9 too;
#    a k-digit decimal number has at most k hex digits only below 10^7 (16^6 > 10^7 > ...): enumerate N < 10^8
found = []
for N in range(10, 10**8):
    h = hex(N)[2:]
    if any(ch > '9' for ch in h): continue
    if str(N)[::-1].lstrip('0') == h.lstrip('0') and is_reversal(N, 10, 16):
        found.append(N)
print("all decimal/hex reversals below 1e8:", found)
assert all(N % 53 == 0 for N in found)
print("every one is divisible by 53")
