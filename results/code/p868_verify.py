"""Puzzle 868: check the five 3x3 magic squares of 101-digit palprimes in results/data/868_squares.txt.
Needs: pip install cypari2 python-flint (both give primality proofs)."""
import cypari2, flint
pari = cypari2.Pari()

blocks = open('results/data/868_squares.txt').read().split('\n\n')
squares = []
for b in blocks:
    nums = [int(x) for x in b.split() if x.isdigit()]
    if nums:
        assert len(nums) == 9
        squares.append([nums[0:3], nums[3:6], nums[6:9]])
assert len(squares) == 5

for sq in squares:
    flat = [x for r in sq for x in r]
    lines = sq + [list(c) for c in zip(*sq)] + [[sq[i][i] for i in range(3)], [sq[i][2 - i] for i in range(3)]]
    assert len({sum(l) for l in lines}) == 1                 # magic: rows, columns, diagonals
    assert len(set(flat)) == 9                               # nine different numbers
    assert all(len(str(x)) == 101 and str(x) == str(x)[::-1] for x in flat)   # 101-digit palindromes
    assert all(pari.isprime(x) for x in flat)                # PARI: proof (APR-CL / ECPP)
    assert all(flint.fmpz(x).is_prime() for x in flat)       # FLINT: independent proof
    print("magic sum", sum(sq[0]))
print("all checks pass")
