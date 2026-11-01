# CLAUDE.md

Goal: beat records and answer open questions on primepuzzles.net. Use reasoning, not brute force.

Key files:
- `README.md`: the results tables and the status of every page.
- `IDEAS.txt`: tricks and open targets. Read it before you pick a puzzle. Add what you learn.
- `FINDINGS.txt`: surprising results. Add new ones.
- `FACTS.md`: facts about primes, each with a formula and a status tag (theorem, checked, computed, literature, conjecture). Add a fact only with a source or a check in `writeup/code/facts_check.py`.

## Language

Write in simple English. This holds for this file, for write-ups, for code comments and for answers in chat.
- Use short sentences and common words. Put one idea in each sentence.
- Use plain words, not jargon. If you need a technical term, explain it the first time you use it.
- Prefer lists and tables to long paragraphs.
- Write facts, not stories. Leave out praise, filler and remarks about yourself.

## Where to find things

- `data/INDEX.txt`: one line per page: `puzzle-NNNN <url> <title>`. Lines with `->` list the attached documents.
- `data/puzzles/puzzle-NNNN.txt` (also `problem-`, `conjecture-`): the full text of a page. NNNN has 4 digits. The page text starts after the `=====` line. Raw HTML is in `data/puzzles/html/`.
- The first line of a page file is the site's HTML title. On some pages it is a stale copy of the previous title (puzzles 314, 319, 344, 355, 393, 543, 768, 821, 849, 900, 986, 1056, 1130, 1233 and conjecture 101). Trust the heading inside the page text.
- `references/*.pdf`, `references/*.txt`: documents linked from the pages, with their text. They include catalogues and tables (for example Keller's `Catalogue.txt` of Brier numbers).
- `references/web/*.txt`: copies of external web pages. Line 1 gives the URL. To find one: `grep -l "^Source: <url>" references/web/*`.
- `references/oeis/`: OEIS data (CC BY-SA 4.0, see its README.md).
  - `puzzle-sequences.txt`: every sequence cited on the pages; `puzzle-sequences-index.txt` says which page cites it.
  - `names.gz`, `stripped.gz`: names and terms of all ~400k sequences. `oeisdata/`: full entries (sparse clone, `clone_oeisdata.sh`). These are git-ignored: rebuild them before you rely on them.
  - Usage: `lookup.py 1,3,13,63` (terms), `lookup.py -w <words>` (names), `lookup.py A001850 --full` (full entry).
- `writeup/NNN.md` (no leading zeros), `writeup/conj-NNN.md`, `writeup/prob-NNN.md`: our write-ups. `writeup/code/` holds the code, `writeup/data/` the data.

## How to attack a puzzle

1. Read the whole page, to the last line. Later contributions often hold the real record. Write down the record and who holds it.
2. Find every page that links to the puzzle, not only the previous and next page: `grep -l -E "[Pp]uzzles? *$N\b|puzz_0*$N\.htm" data/puzzles/*.txt`. Follow-up pages often hold better results.
3. Look for related work: a solved sibling in `writeup/`, `IDEAS.txt`, `FACTS.md`, the OEIS (for any sequence), and the catalogues in `references/`. Puzzles come in families, and a trick that worked once often works again.
4. Explain why the current record works before you write code.
5. Break the problem down to its smallest parts. Check each part on real data before you build on it.
6. Look for structure:
   - Factor the record numbers and check them modulo small numbers (3, 4, 5, 8, 40).
   - Ask whether the property depends only on a remainder. If it does, search only inside that remainder class.
   - Try known record objects first (Prime Pages records, long prime progressions, record polynomials).
7. Try at least three approaches. Write down why each one failed.
8. Use 2 or 3 helper agents at most. Check everything they report yourself.

## General insights

- Many "best found" records are exact optimisation problems: covers, budgets, covering systems, sums of fractions. First shrink the problem by reasoning. Then solve it exactly with an integer program (`scipy.optimize.milp`), a SAT solver (`pysat`) or OR-Tools CP-SAT.
- Ways to shrink a problem:
  - Keep only the largest useful objects. A cover can always use a bigger set instead of a smaller one.
  - Group objects by what they have in common, for example the digit counts of a number.
  - Split the problem by parity or by a remainder.
  - Let objects share their common parts, for example words with the same prefix.
- An LP bound can prove optimality. Give the dual weights as exact fractions and check them in a script. Handle the infinitely many large cases with one extra argument. When the program has few constraints, the best whole-number answer is often the LP bound rounded down.
- SAT solvers are good at word, string and grid-reading conditions. They are bad at "cover everything with at most N sets". For covers, use branch-and-bound with LP bounds.
- For arrangement puzzles (grids, rings), place the hard items first by reasoning and fix them. Then let a local search place the rest. A local search alone cannot build groups of several cells that belong together.
- Check every object on a page against the stated rules with a script. Pages and catalogues contain mistakes.

## Code

- No open-ended hunts. A bounded search is fine once reasoning has made the space small. State every range in the write-up.
- Remove candidates with small factors first (a `numpy` sieve). Test only the survivors with `gmpy2.is_prime`.
- Vectorise with `numpy`. If Python is too slow, write the check in C. Check the machine first (cores, memory, compiler).
- Run every job with `timeout` and `nice`. Do not start long background jobs.
- When you wait for a job, use a check that cannot match itself. `pgrep -f "x"` inside a loop also finds the loop.
- Keep scratch files outside the repo.

## Rules for claims

- Check every result with a second, separate script that recomputes it from scratch (`sympy.isprime`, plain integers).
- Above 2^64, a prime without a primality proof is "PRP". `sympy.isprime` is a proof only below 2^64.
- Results from a solver (integer program, SAT, CP-SAT):
  - The validation script checks every solution.
  - A claim that a value is optimal or impossible counts as "computed" until a second, independent solver or an exact certificate agrees. Name the solver, the time limit and the result in the write-up.
  - Results of a local search or other heuristic are "best found", never "optimal".
- Before you call a result "first" or "new", search the linked pages, the follow-up pages and the catalogues in `references/`.
- Use "proved", "minimal" and "first" only with evidence. Say what you checked and what you assumed.
- Work out the math in full: definitions, every formula derived step by step, every constant explained. Label each step proved, computed, heuristic or empirical, and check each formula on real data.
- Every number in a write-up must come from the validation script or from a derivation in The math.

## Write-up

- Copy the layout of `writeup/1230.md`:
  - **In short:** facts only, no names. What was claimed, asked or held as the record; the proof or construction in a few numbered steps; what it proves or beats. Keep it short.
  - Then: Puzzle link, Status, The question, Record vs ours, How we found it, The math, Validation (a script that checks every claim).
- Link every code or data file you mention, with a path relative to the write-up (`code/x.py`, `data/x.txt`).
- Status line:
  - A beaten record or an answer: add a row to `README.md` under the newest Results heading and update the count in the table heading.
  - A partial result: start with `open`. No row in `README.md`.
  - The page already has a proved final answer: start with `closed` and cite the proof.
  - A conjecture: start with proved, disproved, explained or answered.
- Then run `python3 writeup/code/make_status.py` and `python3 writeup/code/stamp.py`.
- Make one commit per puzzle: `Add puzzle NNN: ...` (or `Add conjecture NNN: ...`).

## Environment

The scripts need `sympy`, `gmpy2` and `numpy`; some also need `scipy`, `python-sat`, `ortools` or `cypari2`. Use a local virtualenv (git ignores `.venv/`):

```
python3 -m venv .venv && .venv/bin/pip install sympy gmpy2 numpy scipy python-sat ortools
```
