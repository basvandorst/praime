#!/bin/sh
# Recreate references/oeis/oeisdata: a shallow, sparse clone of github.com/oeis/oeisdata (git-ignored).
# Only the .seq files of the A-numbers cited in data/puzzles are checked out (about 80 MB instead of gigabytes).
# Add more:  cd references/oeis/oeisdata && echo /seq/A000/A000040.seq >> .git/info/sparse-checkout && git checkout main
set -e
cd "$(dirname "$0")"
[ -d oeisdata ] || git clone -q --depth 1 --filter=blob:none --no-checkout https://github.com/oeis/oeisdata.git oeisdata
cd oeisdata
git sparse-checkout init --no-cone
{ echo /LICENSE; echo /README.md
  grep -ohE "A[0-9]{6}" ../../../data/puzzles/*.txt | sort -u | awk '{printf "/seq/%s/%s.seq\n", substr($0,1,4), $0}'; } > .git/info/sparse-checkout
git checkout -q main
echo "$(find seq -name '*.seq' | wc -l) sequences checked out"
