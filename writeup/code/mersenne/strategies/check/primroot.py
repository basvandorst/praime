"""Independent check of the primroot note: least primitive roots of the odd primes up to 500,000 with sympy."""
import sympy, collections
c = collections.Counter(); mx = (0, 0); sq = 0; m8 = 0
for p in sympy.primerange(3, 500001):
    g = sympy.primitive_root(p, smallest=True)
    c[g] += 1
    if g > mx[0]: mx = (g, p)
    if int(g ** .5) ** 2 == g: sq += 1
    if g == 2 and p % 8 in (1, 7): m8 += 1
tot = sum(c.values())
print('odd primes', tot, '; g=2:', c[2], f'{100*c[2]/tot:.2f}%', '; g=3:', f'{100*c[3]/tot:.2f}%', '; max', mx, '; square g:', sq, '; g=2 with p = +-1 mod 8:', m8)
