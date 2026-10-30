# CLAUDE.md

Goal: beat records and answer open questions on primepuzzles.net. Use reasoning, not brute force. `README.md` has the results tables.

## Where to find things

- `data/INDEX.txt`: one line per page: `puzzle-NNNN <url> <title>`. Lines with `->` give the attached documents.
- `data/puzzles/puzzle-NNNN.txt` (also `problem-`, `conjecture-`): full page text. The number has 4 digits. The page text starts after the `=====` line. Raw HTML is in `data/puzzles/html/`.
- `references/*.pdf` and `references/*.txt`: documents linked from the pages, with extracted text.
- `references/web/*.txt`: copies of external web pages. Line 1 gives the URL. To find one: `grep -l "^Source: <url>" references/web/*`.
- `results/NNN.md`: our write-up per puzzle (no leading zeros). `results/code/` and `results/data/` hold large code and data.

Some titles do not match the file number (`puzzle-0355.txt` is "Puzzle 354"). Check the title.

## How to crack a puzzle

1. Read the full page, to the last line. Later contributions often hold the real record. Write down the record and who holds it.
2. Find a solved sibling in `results/`. Puzzles come in families: prime-generating formulas, magic squares, digit rules, σ/φ equations, covering sets. Use the trick that worked there.
3. Explain why the current record works before you write code.
4. Factor the record numbers. Check them mod 3, 4, 5, 8 and 40. Look for a pattern.
5. Turn the condition into a remainder. If the property depends only on p mod m, build a formula that always gives that remainder. Then search only that class.
6. Try known record objects first: Prime Pages records, AP-27, Dress–Landreau polynomials.
7. Try at least three approaches. Write down why each one failed.
8. Search only in bounded ranges. Sieve small factors first. Then test the survivors with `gmpy2.is_prime`. Use numpy or C for hot loops.
9. Check every result again with a second, independent script.

## Rules

- Limit every job with `timeout` and `nice`. Do not start long background jobs.
- Above 2^64, a prime without a primality proof is "PRP". (`sympy.isprime` is a proof only below 2^64.)
- Do not overclaim. Use "proved", "minimal" and "first" only with evidence. Say what you checked and what you assumed.
- Put scratch files outside the repo.

## Write-up

Copy the layout of `results/1230.md`: In short, Puzzle link, Status, The question, Record vs ours, How we found it, The math, Verify it yourself (a script that checks every claim). Add a row to `README.md` and update the counts in its first lines. Make one commit per puzzle: `Add puzzle NNN: ...`.

## Environment

The scripts need `sympy`, `gmpy2` and `numpy`. If they are missing, make a local virtualenv (git ignores `.venv/`):

```
python3 -m venv .venv && .venv/bin/pip install sympy gmpy2 numpy
```

Check what the machine has (cores, memory, a C compiler) before you plan a search. Use C only if a compiler is installed.
