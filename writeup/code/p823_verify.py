# Puzzle 823: K(S) = number of distinct primes among the substrings of S (digits 1, 3, 7, 9), read left to right.
# Independent evaluation (sympy below 2^64, PARI isprime above) of the strings on the page and of ours.
# Run from the repo root: .venv/bin/python3 writeup/code/p823_verify.py
from sympy import isprime
import cypari2
pari = cypari2.Pari()
def K(s):
    found = set()
    for i in range(len(s)):
        for j in range(i + 1, len(s) + 1):
            t = s[i:j]
            if t in found: continue
            v = int(t)
            if (isprime(v) if v < 2**64 else pari.isprime(v) == 1): found.add(t)
    return len(found)
page = [("Rivera L=20", "19717393139337391999", 68), ("Vantieghem L=20", "91131179719337391999", 73), ("Kamenetsky L=20", "39311171399719337397", 74),
        ("Vantieghem L=40", "7139971933199733313711131179719337391999", 156), ("Kamenetsky L=40", "3919973911171131733313997933739717193337", 166),
        ("Hurter L=40", "1733391137971933199739979337979931793719", 183),
        ("Kamenetsky L=80", "19133139391911371797171391397919391117191773719393977397333113933733937199937799", 354),
        ("Hurter L=80", "13999917333119797193319973711793793333133333373919393911317773131939913739397117", 431)]
ours = [("ours L=20 (annealing)", "93719333139979337397", 74), ("ours L=40 (annealing, random start)", "9917331113193937393199733313939191719333", 162)]
for name, s, k in page + ours:
    assert set(s) <= set("1379") and len(s) in (20, 40, 80)
    got = K(s); print(f"{name}: K = {got} (claimed {k})"); assert got == k, name
print("all checks passed")
