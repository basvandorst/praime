# A121805: a(1) = 1, a(n+1) = smallest number > a(n) with a(n+1) - a(n) = 10*(last digit of a(n)) + (first digit of a(n+1)).
from collections import Counter
def first(x): return int(str(x)[0])
a = 1; seq = [1]
while True:
    L = a % 10; nxt = None
    for F in range(1, 10):
        c = a + 10 * L + F
        if first(c) == F: nxt = c; break
    if nxt is None: break
    a = nxt; seq.append(a)
print("lengte", len(seq), "laatste", seq[-1])
D = Counter(b - a for a, b in zip(seq, seq[1:]))
print("D=5,55:", D[5], D[55], " D=15,65:", D[15], D[65], " D=4,24..84:", [D[x] for x in (4, 24, 44, 64, 84)])
# block entries: first term of each maximal block with constant (number of digits, first digit)
entries = []
for i in range(1, len(seq)):
    if (len(str(seq[i])), first(seq[i])) != (len(str(seq[i-1])), first(seq[i-1])):
        entries.append((len(str(seq[i])), first(seq[i]), seq[i - 1], seq[i], seq[i] % 10))
for F in (4, 5):
    print(f"blokken met eerste cijfer {F}: (cijfers, vorige term, eerste term, laatste cijfer van eerste term)")
    for e in entries:
        if e[1] == F: print("   ", e[0], e[2], e[3], e[4])

# predict every count from the block structure: a block (digits j, first digit F) entered with last digit L0
# uses the D values 10*L + F for L in the orbit L0, L0+F, L0+2F, ... (mod 10); its length ~ size / mean step,
# shared equally by those D values.
pred = Counter()
bounds = entries + [(None, None, seq[-1], None, None)]
for i, (j, F, prev, firstterm, L0) in enumerate(entries):
    lo = firstterm; hi = min(bounds[i + 1][3] if bounds[i + 1][3] else seq[-1], (F + 1) * 10**(j - 1))
    orbit = []; L = L0
    while L not in orbit: orbit.append(L); L = (L + F) % 10
    Ds = [10 * l + F for l in orbit]; mean = sum(Ds) / len(Ds)
    for d in Ds: pred[d] += (hi - lo) / mean / len(Ds)
worst = max(abs(pred[d] - D[d]) / D[d] for d in D if D[d] > 1000)
print("voorspeld vs echt:", {d: (round(pred[d]), D[d]) for d in (5, 55, 2, 8, 14, 1, 45, 52, 35, 6)})
print(f"grootste relatieve fout voor D met >1000 voorkomens: {worst:.4f}")
print("D die nooit voorspeld worden:", sorted(d for d in range(100) if pred[d] < 0.5 and d % 10))
