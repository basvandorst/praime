#!/bin/bash
# The later runs of machine 23, as they were made on 18 and 19 November 2026 (each line one run; not meant to be run as a whole).
PY=/home/bas/praime/.venv/bin/python3
# random primes of sizes that the rule was not made from (seed 20261118 + digits + 1000 K, K the third argument)
$PY inherit_big.py 65 60   | grep digits > inherit_big_65.out
$PY inherit_big.py 70 40   | grep digits > inherit_big_70.out
$PY inherit_big.py 75 40   | grep digits > inherit_big_75.out
$PY inherit_big.py 75 40 1 | grep digits > inherit_big_75b.out   # drawn because the first sample of 75 digits stood out
$PY inherit_big.py 80 24   | grep digits > inherit_big_80.out
$PY inherit_big.py 85 20   | grep digits > inherit_big_85.out
for spec in "11 3000" "13 3000" "14 3000" "16 3000" "17 3000" "18 3000" "19 3000" "22 3000" "28 2000" "32 1000" "45 300"; do set -- $spec; $PY inherit_big.py $1 $2 | grep digits; done > inherit_big_more.out
# the scans
./inherit 4200000000 > inherit_c_4e9.out ; ./inherit 4200000000 0 5 inherit_deep_4e9.txt ; ./inherit 1000000000 1 > inherit_c_minus_1e9.out
./inherit_second 2100000001 4200000000 > inherit_second_4e9.out
./inherit_second 4200000001 8400000000 > inherit_second_8e9.out      # then the lines "deep 5 ..." moved to inherit_deep_8e9.txt
./inherit_third 8400000001 16800000000 > inherit_third_17e9.out      # then the lines "deep 5 ..." moved to inherit_deep_17e9.txt
# the lemma instead of a scan
$PY inherit_lemma.py > inherit_lemma.txt ; $PY inherit_lemma.py 16800000000 inherit_deep_8e9.txt > inherit_lemma_2.txt
$PY inherit_lemma_chain.py 4200000000 inherit_deep_4e9.txt > inherit_lemma_chain_4e9.txt
$PY inherit_lemma_chain.py 8400000000 inherit_deep_8e9.txt > inherit_lemma_chain_8e9.txt
$PY inherit_lemma_chain.py 16800000000 inherit_deep_17e9.txt > inherit_lemma_chain_17e9.txt
# the search by construction for primes of 7 rounds
$PY inherit_deep.py 2000 200000 > inherit_deep_first.out             # stopped by hand in stage 2; stage 1 wrote inherit_deep_six.txt
$PY inherit_deep.py stage2 320000000000000 > inherit_deep2.out
$PY inherit_deep.py 100000000 0 39000000000000 > inherit_deep_wide.out   # stage 1 with every even s up to min(r, 3.9e13 / r); writes inherit_deep_six_wide.txt
$PY inherit_deep.py stage2 77000000000000 inherit_deep_six_wide.txt > inherit_deep_wide2.out
# the rounds with K sources of helpers
for k in 1 2 3 5; do $PY inherit_sources.py $k | grep "^{" > inherit_sources_$k.out; done
