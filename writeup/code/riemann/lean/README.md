# The Lean proof of RH and GRH for Dirichlet L-functions (November 2026)

Static checks of the announced Lean proof, the page that shows them (`prime38.html`), and the prime sum of the proof (FACTS.md, facts 2047 and 2048). The archive itself (`lean.tar.gz.part01` to `part14` from https://alpo.ge, 239,275,554 bytes, SHA-256 `9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc`) is not kept here. Every script streams it from a folder given on the command line.

## The archive

| file | what it does |
|---|---|
| `scan_lean.py` | streams every `.lean` file and counts `axiom`, `sorry`, `native_decide`, `implemented_by`, `unsafe`, environment edits, macros and `set_option` (comments removed) |
| `graph_scan.py` | one pass: per file the bytes, lines, declarations and imports → `graph.json` (88 MB, not kept) |
| `graph_analyse.py` | the import cone of `Main` (all 175,187 files), no cycles, no missing imports, the longest chain (202 steps) |
| `page_data.py` | the small data block of the page → `data38.json` |
| `ref38.py` | an independent recount (other reader, `splitlines`, token parser, BFS, Kahn layers); needs `lean_list.txt` = `cat parts \| tar -tvz` → `ref38.json` |
| `compare38.py` → `compare38.txt` | page data against `ref38.json` and the page's zeros and ζ values against mpmath: 0 mismatches |
| `verify_review.py` | recount of the points of the independent review: home-made kernel tactics (717 uses in 86 files), literature `Prop`s (2,020 without parameters), `Counterexamples` folders |
| `graph_scan_cat.py` | the multi-part reader shared by the scans |

## The page

`prime38.src.html` + `data38.json` → `gen38.py` → `prime38.html`. `t38.js` runs it in headless Chrome. Sections: ζ in the strip and Z(t); Weil's positivity from the primes (`weil_ref.py` → `weil_ref.json`, `weil_compare.txt`: the page's terms agree with mpmath to 5·10⁻¹²); the primes in a tube; the treemap and the import layers; the last Lean steps; who checked what.

## The tube (fact 2048)

| file | what it does |
|---|---|
| `tube.py N` | G(y) = Σ Λ(n)(y/n − 1) − y ln y + (1 + γ)y at all prime powers ≤ N and at the maxima in between; `tube_1e7.txt`, `tube_1e8.txt` |
| `tube_zeros.py` | G(y) from the primes against ln 2π − T(y) + Σ over 400 zeros (`weil_ref.json`): the formula of fact 2048 |
| `check2048.py` | the assertions of fact 2048 (also in `writeup/code/facts_check.py`) |
| `tube_compare.txt` | the page's tube values against these scripts |

Run with the repository's virtualenv, for example `(ulimit -v 2500000; timeout 900 nice ../../../../.venv/bin/python -I tube.py 1e8)`.

## The short edition (github.com/lsernau/grh-lean, 25 November 2026)

A second Lean project proves the same five theorems in about 300,000 lines (readable edition, `grh/`) or 76,000 lines (`grh-compact/`), on Lean v4.35.0, with a reading guide (`grh/GUIDE.md`), an index of results with sources (`INDEX.md`) and a list of references. Its README says: based on Alpöge's proof, many steps by other arguments; machine checks only, no mathematician has read it through. Read here (commit 72fead10, cloned, nothing built): README, GUIDE, the first tables of INDEX, and the heads of `New/Reduction/{MainInequality, PrimeSums, OneSidedBound}`, `New/Reduction/Zeros/{Criterion, PrimeError}`, `Lib/Analytic/QuadraticTwist`, `New/Derivatives/Divisibility/Count`.

What the code says about primes:
- **The whole proof ends in one inequality on primes.** For q₀ a multiple of 18480 and a prime d ≡ q₀/2 − 1 mod q₀: `R_q(x) = Σ_(n≤x, n≡1 mod q) Λ(n)·(1 − (n/d))·(1/n − 1/x) − ln x/φ(q) − c_q ≤ (64φ(q₀) + 3)·x^(−1/2+6ε)` for large integers x (`eventually_primeRemainder_natCast_le`). One side only.
- **The criterion** (`LFunction_ne_zero_of_primeRemainder_le`, 91 lines): if R_q(x) ≤ x^(β−1) eventually for every prime d of one class mod q₀ (8 | q₀, d ≡ 7 mod 8), then no L(s, ψ), ψ mod q₀, vanishes in Re s > β. The weight 1 − (n/d) makes the zeros of L(s, ψ) enter with a plus sign and those of the twists L(s, ψ·(·/d)) with a minus sign; a prime d whose twists have no zero in a given box exists by a zero-density count that needs only Pólya–Vinogradov (`Lib/Analytic/QuadraticTwist`); then Landau's lemma on a transform with non-negative residues.
- **Where the prime sum comes from:** size against divisibility. A Gram determinant of theta lifts is an algebraic number, so its sizes at all places add up to at least 0 (product formula, `TwistedLattice.covolume_add_index_nonneg`). The size at infinity of a derivative of order j is ln j! + j·ln(M/Ω); the divisibility at p of the same derivatives is a Legendre-type count; their difference is the prime sum Σ Λ(n)(1 − (n/d))·(y/n − κ)₊ less its expected value. The main amounts cancel and the remainder is bounded above. The literature the code names for this idea: Gelfond–Schnirelman (Pritsker 2005), Granville–Stark, Bombieri's proof for curves.
- **The exponent 1/2** comes from the number of lifts m, between √x and 2√x: the error per lift is log⁶ b / m.
- **The constant** Ω = exp(γ − Re L′/L(1, χ_K))·2π/rd_E is a Chowla–Selberg period; c₁ = Re G₁(1) cancels between the period and the expected value of the prime sums.

Own check on primes ([rq_tube.py](rq_tube.py) → [rq_tube_out.txt](rq_tube_out.txt)), for the smallest pair the criterion allows, q = 8 and d = 7 (the proof itself needs 18480 | q₀): c₈ = −0.6944767878 from L′/L(1, χ_D), D = −4, 8, −8, −7, 28, −56, 56; with GRH for Dirichlet L-functions (fact 2047) the exact formula gives |R₈(x)|·√x ≤ C = 0.9335886 (the sum of 1/(ρ(1 − ρ)) over the zeros of the eight L-functions, divided by 4). On the 720,656 prime powers n ≡ 1 mod 8 with (n/7) ≠ 1 up to 10⁸, R₈(x)·√x stays between −0.377 and +0.419 (0.45 of C), and at x = 10⁴, 10⁶, 10⁸ the primes give +0.2605, −0.2113, −0.2005 against +0.2711, −0.2091, −0.2002 from the zeros below 600. So the two-sided truth is a band of fixed width in units of 1/√x, where the proof shows one side with x^(6ε) to spare.
