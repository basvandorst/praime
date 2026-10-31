# Local OEIS reference

Content from The On-Line Encyclopedia of Integer Sequences (OEIS), OEIS Foundation Inc., https://oeis.org.
Licence: Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0),
https://creativecommons.org/licenses/by-sa/4.0/ (see also https://oeis.org/LICENSE).
The OEIS LICENSE file says: "The OEIS content is copyrighted by the OEIS Foundation, Inc.
The content is made available under the terms of the CC-BY-SA-4.0 license".
The two text files here are excerpts, modified (trimmed), retrieved 2026-10-31, and are shared under the same licence.

## Files

- `puzzle-sequences.txt`: the 613 OEIS sequences cited in `data/puzzles/*.txt`. One block per sequence, sorted by A-number:
  A-number, name, the whole term list (%S/%T/%U) on one line, `F:` lines (all %F formulas), `K:` keywords,
  `offset:`, and at most 3 `C:` comments cut to 300 characters. Programs, references, links, examples and
  cross-references are dropped.
- `puzzle-sequences-index.txt`: one line per A-number: `A-number<TAB>name<TAB>puzzle files citing it`.
- `lookup.py`: offline search. `python3 lookup.py 2,3,5,7,11,13` lists the sequences that contain these
  terms consecutively; `python3 lookup.py A005384` prints a name.
- `stripped.gz`, `names.gz` (not in git): the OEIS daily bulk files, https://oeis.org/stripped.gz and
  https://oeis.org/names.gz. `lookup.py` needs them and says how to download them when they are missing.
- `oeisdata/` (not in git): sparse, shallow clone of https://github.com/oeis/oeisdata with the `.seq` files
  (OEIS internal format) of the cited sequences. `clone_oeisdata.sh` recreates it.

## How the text files were made

1. `grep -ohE "A[0-9]{6}" data/puzzles/*.txt | sort -u` gives 649 hits.
2. 25 hits have no OEIS entry (A000000, A999999, ...): they are digit strings in puzzle text.
   11 more are real A-numbers that occur only as the head of a longer digit string (such as A0000058
   in problem-0083), so they are not citations. Both groups are left out; the header of
   `puzzle-sequences.txt` lists them.
3. The remaining 613 entries were read from `oeisdata/seq/AXXX/AXXXXXX.seq`
   (oeisdata commit c24b7a0, Oct 31 2026) and trimmed as described above. The files in the index come from
   `grep -lE "ANNNNNN([^0-9]|$)" data/puzzles/*.txt`.

## Everyday use

```
python3 references/oeis/lookup.py 1,3,13,63,321          # which sequences contain these terms (all ~400k, offline)
python3 references/oeis/lookup.py -w sophie germain      # names containing these words
python3 references/oeis/lookup.py A001850 --full         # whole entry with formulas (fetched into the clone if missing)
python3 references/oeis/lookup.py --update               # refresh the clone; also tries the two .gz files
```

oeis.org sometimes answers scripts with a Cloudflare bot challenge (HTTP 403). Then `--update` keeps the old
`.gz` files and says so; download them in a browser instead. The clone updates from GitHub and is not affected.
Rebuild the excerpt with `CAP=1000 python3 references/oeis/build_extract.py` (it reads `notice.txt`).
