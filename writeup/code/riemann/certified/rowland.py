import math
a, n, firsts = 7, 1, {}
while len(firsts) < 12 and n < 10 ** 7:
    n += 1; g = math.gcd(n, a); a += g
    if g > 1 and g not in firsts: firsts[g] = n
print('Rowland: first appearance of each new prime difference (prime: step): ' + ', '.join(f'{p}: {s}' for p, s in sorted(firsts.items(), key=lambda t: t[1])))
