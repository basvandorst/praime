# L(3*2^n + 1) vs L(2^(n+2) + 1): equal unless the run reaches a k whose own path needs e(k) >= n halvings.
def H(x):
    s = 0
    while x != 1:
        x = x // 2 if x % 2 == 0 else 3 * x + 1; s += 1
    return s
def eo(k):
    e = o = 0
    while k != 1:
        if k % 2: k = 3 * k + 1; o += 1
        else: k //= 2; e += 1
    return e, o
def L(x, cap=10**5):
    h = H(x); r = 1
    while r < cap and H(x + r) == h: r += 1
    return r
exc, pred = [], []
for n in range(3, 141):
    a, b = L(3 * 2**n + 1), L(2**(n + 2) + 1)
    if a != b: exc.append(n)
    # predicted: identity can fail only if some k <= max(a, b) + 1 has e(k) >= n
    if max(eo(k)[0] for k in range(1, max(a, b) + 2)) >= n: pred.append(n)
print("n in 3..140 with L(3*2^n+1) != L(2^(n+2)+1):", exc)
print("n where the run reaches a k with e(k) >= n:   ", pred)
print("e(27), o(27) =", eo(27))
assert set(exc) <= set(pred)
