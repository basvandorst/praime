# CLAUDE.md

Goal: beat records and answer open questions on primepuzzles.net. Use reasoning, not brute force. `README.md` has the results tables. `FINDINGS.txt` lists surprising results. `IDEAS.txt` lists tricks and open targets: read it before you pick a puzzle, and update it when you learn something new. `FACTS.md` lists facts about primes with formulas and status tags (theorem, checked, computed, literature, conjecture); add a fact there only with a source or a check in `results/code/facts_check.py`.

## Where to find things

- `data/INDEX.txt`: one line per page: `puzzle-NNNN <url> <title>`. Lines with `->` give the attached documents.
- `data/puzzles/puzzle-NNNN.txt` (also `problem-`, `conjecture-`): full page text. The number has 4 digits. The page text starts after the `=====` line. Raw HTML is in `data/puzzles/html/`.
- `references/*.pdf` and `references/*.txt`: documents linked from the pages, with extracted text.
- `references/web/*.txt`: copies of external web pages. Line 1 gives the URL. To find one: `grep -l "^Source: <url>" references/web/*`.
- `results/NNN.md`: our write-up per puzzle (no leading zeros). `results/code/` and `results/data/` hold large code and data.
- `FACTS.md`: the reference list of facts about primes, each with its formula and how it is known: residue rules (section 1), counts (2), gaps (3), patterns and runs of consecutive primes (4), prime-producing formulas (5), record lengths (8), facts from our write-ups (9) and from the archive (10), open problems (11). Read the relevant section before you attack a puzzle.
- `writeup/NNN.md`: attempts on puzzles that are still open. `writeup/INDEX.md` lists all puzzles with their status (open, beaten, answered).

Some titles do not match the file number (`puzzle-0355.txt` is "Puzzle 354"). Check the title.

## How to crack a puzzle

- Read the full page, to the last line. Later contributions often hold the real record. Write down the record and who holds it. An audit withdrew seven false claims that came from reading only part of a page (855, 872, 933, 969, 854, 1136, and k = 14 in 522).
- Find a solved sibling in `results/` and check `IDEAS.txt` and `FACTS.md`. Puzzles come in families. A trick that worked once often works again.
- Explain why the current record works before you write code. Then propose a construction.
- Factor the record numbers with `sympy.factorint`. Check them mod 3, 4, 5, 8 and 40. Patterns show up fast.
- Turn the condition into a remainder. Ask whether the property ("even period", "unlucky", "no digit 6") depends only on p mod m. If it does, build a formula that always gives that remainder, and search only inside that class.
- Try known record objects first: Prime Pages records, PrimeGrid's AP-27, Dress–Landreau polynomials (`FACTS.md` sections 8–10 list them). They often win outright.
- Try at least three approaches. Write down why each one failed.
- Split work only within limits: 2 or 3 agents at most, for example one per page range for scouting or one per hard puzzle. Re-check everything they report. One agent invented an explanation.

## Code

- No open-ended hunts. A bounded exhaustive search is fine once reasoning has made the space small (81, 106, 75 rest on one). State every range in the write-up.
- Sieve before testing. Remove candidates with small factors with a `numpy` boolean sieve, then run `gmpy2.is_prime` only on the survivors.
- Vectorise exhaustive checks with `numpy`. If Python is still too slow, write the check in C.
- Limit every job with `timeout` and `nice`. Do not start long background jobs.
- Verify with a second, separate script. Recompute the result from scratch with `sympy.isprime` and plain integers before you call it a record. Several "records" did not survive this.

## Rules

- Above 2^64, a prime without a primality proof is "PRP". (`sympy.isprime` is a proof only below 2^64.) Puzzles 419, 793 and 908 had to be relabelled for this.
- Do not overclaim. Use "proved", "minimal" and "first" only with evidence. Say what you checked and what you assumed.
- Keep a running log of checks. Put scratch files outside the repo.
- Add what you learn to `IDEAS.txt` and `FINDINGS.txt`.

## Write-up

Copy the layout of `results/1230.md`: In short, Puzzle link, Status, The question, Record vs ours, How we found it, The math, Verify it yourself (a script that checks every claim). Add a row to `README.md` and update the counts in its first lines. An attempt that does not beat or answer the puzzle goes to `writeup/NNN.md` instead. After either, run `python3 writeup/make_index.py`. Make one commit per puzzle: `Add puzzle NNN: ...`.

## Environment

The scripts need `sympy`, `gmpy2` and `numpy`. If they are missing, make a local virtualenv (git ignores `.venv/`):

```
python3 -m venv .venv && .venv/bin/pip install sympy gmpy2 numpy
```

Check what the machine has (cores, memory, a C compiler) before you plan a search. Use C only if a compiler is installed.
