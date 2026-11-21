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
