# CLAUDE.md

Goal: beat records and answer open questions on primepuzzles.net. Use reasoning, not brute force. `README.md` has the results tables. `FINDINGS.txt` lists surprising results. `IDEAS.txt` lists tricks and open targets: read it before you pick a puzzle, and update it when you learn something new. `FACTS.md` lists facts about primes with formulas and status tags (theorem, checked, computed, literature, conjecture); add a fact there only with a source or a check in `writeup/code/facts_check.py`.

## Where to find things

- `data/INDEX.txt`: one line per page: `puzzle-NNNN <url> <title>`. Lines with `->` give the attached documents.
- `data/puzzles/puzzle-NNNN.txt` (also `problem-`, `conjecture-`): full page text. The number has 4 digits. The page text starts after the `=====` line. Raw HTML is in `data/puzzles/html/`.
- `references/*.pdf` and `references/*.txt`: documents linked from the pages, with extracted text.
- `references/web/*.txt`: copies of external web pages. Line 1 gives the URL. To find one: `grep -l "^Source: <url>" references/web/*`.
- `references/oeis/`: OEIS data (CC BY-SA 4.0, see its README.md).
  - `puzzle-sequences.txt` holds name, terms and formulas of every sequence cited on the puzzle pages, and `puzzle-sequences-index.txt` says which pages cite it.
  - `oeisdata/` is a sparse clone of github.com/oeis/oeisdata with the full `.seq` entries, made with `clone_oeisdata.sh`; add a sequence via `.git/info/sparse-checkout`.
  - `names.gz` and `stripped.gz` hold the names and terms of all ~400k sequences; `lookup.py` searches them offline.
  - The clone and the `.gz` files are git-ignored. Rebuild them before relying on them.
  - Usage: `lookup.py 1,3,13,63` (terms), `lookup.py -w <words>` (names), `lookup.py A001850 --full` (entry with formulas).
  - For every puzzle that involves a sequence, look it up first: OEIS often has more terms, a formula or a sibling.
- `writeup/NNN.md`: our write-up per puzzle (no leading zeros), whether it beat the record or is still open; `writeup/conj-NNN.md` for a conjecture page, `writeup/prob-NNN.md` for a problem page. The Status line says which. `writeup/code/` holds the code and `writeup/data/` the large data. The Status section of `README.md` lists every puzzle, conjecture and problem as beaten, solved, partially solved or unsolved.
- `FACTS.md`: the reference list of facts about primes, each with its formula and how it is known: residue rules (section 1), counts (2), gaps (3), patterns and runs of consecutive primes (4), prime-producing formulas (5), record lengths (8), facts from our write-ups (9) and from the archive (10), open problems (11). Read the relevant section before you attack a puzzle.

The first line of a page file is the site's HTML title, and on 14 puzzle pages (314, 319, 344, 355, 393, 543, 768, 821, 849, 900, 986, 1056, 1130, 1233) and on conjecture 101 it is a stale copy from the previous page: `puzzle-0355.txt` starts with "Puzzle 354" but is Puzzle 355. Trust the heading inside the page text.

## How to crack a puzzle

- Read the full page, to the last line. Later contributions often hold the real record. Write down the record and who holds it. An audit withdrew seven false claims that came from reading only part of a page (855, 872, 933, 969, 854, 1136, and k = 14 in 522).
- Find a solved sibling in `writeup/` and check `IDEAS.txt` and `FACTS.md`. Puzzles come in families. A trick that worked once often works again.
- Explain why the current record works before you write code. Then propose a construction.
- Factor the record numbers with `sympy.factorint`. Check them mod 3, 4, 5, 8 and 40. Patterns show up fast.
- Turn the condition into a remainder. Ask whether the property ("even period", "unlucky", "no digit 6") depends only on p mod m. If it does, build a formula that always gives that remainder, and search only inside that class.
- Try known record objects first: Prime Pages records, PrimeGrid's AP-27, Dress–Landreau polynomials (`FACTS.md` sections 8–10 list them). They often win outright.
- Try at least three approaches. Write down why each one failed.
- Split work only within limits: 2 or 3 agents at most, for example one per page range for scouting or one per hard puzzle. Re-check everything they report. One agent invented an explanation.

## Insights that keep paying off

- Many "best found" records are exact optimisation problems in disguise (covers, digit budgets, covering systems, Egyptian fractions). Shrink the space by reasoning first, then solve exactly: `scipy.optimize.milp` (HiGHS), `pysat` (CaDiCaL, Glucose), OR-Tools CP-SAT. Puzzles 1276, 963 and 68, and problems 76 and 79, fell this way.
- Shrink before solving. Use maximal objects only (963: maximal progressions), group by invariants (1276: digit-count vectors; 207: digit counts plus order of first appearance), split by parity or residue (68: Brier coverings become two covering problems in t), and share prefixes (896: subsequence SAT over a trie).
- LP duality gives optimality certificates a script can check: weights in exact fractions, plus one argument for the infinite tail (1276: every prime of 8+ digits weighs at least 1). With few constraints the integer optimum is often ⌊LP⌋.
- SAT is strong on word, subsequence and grid-reading constraints (771, 896) and weak on set cover with a cardinality bound. For set cover use branch-and-bound with LP bounds.
- For arrangements with a local cost, place the expensive items by reasoning and freeze them, then let a local search fill in the rest. Plain annealing cannot build motifs of several cells (1277: lonely numbers on the border, each prime surrounded by its multiples).
- Recheck every object on the page against the stated rules with a script. Two of the grids in 771 broke the straight-line rule, and five entries of Keller's catalogue give the wrong P.
- Covering sets: a prime in both sets has even order m and covers n ≡ a and n ≡ a + m/2; with 3 shared, the others need m ≡ 2 (mod 4).

## Code

- No open-ended hunts. A bounded exhaustive search is fine once reasoning has made the space small (81, 106, 75 rest on one). State every range in the write-up.
- Sieve before testing. Remove candidates with small factors with a `numpy` boolean sieve, then run `gmpy2.is_prime` only on the survivors.
- Vectorise exhaustive checks with `numpy`. If Python is still too slow, write the check in C.
- Limit every job with `timeout` and `nice`. Do not start long background jobs.
- Verify with a second, separate script. Recompute the result from scratch with `sympy.isprime` and plain integers before you call it a record. Several "records" did not survive this.
- Solver results (ILP, SAT, CP-SAT): a solution is checked by the verify script; an optimality or impossibility claim counts as "computed" until a second, independent solver or a certificate agrees (an LP dual in exact fractions, a separate branch-and-bound, a second program). Give the solver, the time limit and the status in the write-up. Results of annealing or other heuristics are "best found", never "optimal".
- Before you claim "first" or "new", also search the archived catalogues and tables in `references/` (for example Keller's `Catalogue.txt` for Brier numbers), not only the linked pages.
- Every number in a write-up (counts, bounds, "about X") must come out of the validation script, or be derived step by step in The math.
- Wait for background jobs with a check that cannot match itself: `pgrep -f "x"` inside `until` also finds the loop's own command line.

## Rules

- Above 2^64, a prime without a primality proof is "PRP". (`sympy.isprime` is a proof only below 2^64.) Puzzles 419, 793 and 908 had to be relabelled for this.
- Work out the math in every analysis, down to the smallest detail: formal definitions, every formula derived step by step, every constant explained, each step labelled proved, heuristic or empirical, and each formula checked numerically against real data. This holds for write-ups and for answers in chat alike.
- Do not overclaim. Use "proved", "minimal" and "first" only with evidence. Say what you checked and what you assumed.
- Keep a running log of checks. Put scratch files outside the repo.
- Add what you learn to `IDEAS.txt` and `FINDINGS.txt`.

## Write-up

Link every code or data file a write-up mentions, with a path relative to the write-up (`code/x.py`, `data/x.txt`); paths inside command blocks need a link elsewhere in the same write-up. Copy the layout of `writeup/1230.md`: In short (facts only, no contributor names: what was claimed, asked or held as record; the proof or construction in a few numbered steps; what it proves or beats), Puzzle link, Status, The question, Record vs ours, How we found it, The math, Validation (a script that checks every claim). Add a row to `README.md` under its newest dated Results heading (`## Results DD-MM-YYYY`: the planned release date picked by the user, not the day of the work; conjectures go in a Conjectures table) and update the counts line under that heading. For a conjecture, start the Status line with proved, disproved, explained or answered: the index script reads those first words. An attempt that does not beat or answer the puzzle uses the same layout in `writeup/NNN.md`, with a Status line that starts with open and no row in `README.md`. If the page already holds a proved final answer (proved optimal, unique or complete), start the Status line with closed and cite the proof; the Status section lists it as solved, no improvement possible. After either, run `python3 writeup/code/make_status.py` to rebuild the Status section of `README.md`, then `python3 writeup/code/stamp.py`. Make one commit per puzzle: `Add puzzle NNN: ...` (or `Add conjecture NNN: ...`).

## Environment

The scripts need `sympy`, `gmpy2` and `numpy`. If they are missing, make a local virtualenv (git ignores `.venv/`):

```
python3 -m venv .venv && .venv/bin/pip install sympy gmpy2 numpy
```

Check what the machine has (cores, memory, a C compiler) before you plan a search. Use C only if a compiler is installed.
