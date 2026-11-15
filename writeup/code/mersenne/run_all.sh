#!/bin/bash
# Reruns every experiment of MERSENNE.md, in order. Needs sympy, gmpy2, numpy, scipy and Node 22. About two hours on two cores
# (most of it the P-1 runs in touch/). The .npy files are caches and stay out of git.
set -e
cd "$(dirname "$0")"
P=${PY:-python3}
run() { (ulimit -v 1500000; timeout 4000 nice "$@"); }

# 1. Lucas-Lehmer on every prime below 10^4; hit rate after trial factoring; the heuristic count
run $P ll_check.py 10000 > ll_check.txt
run $P tf.py 10000 100000
run $P tf.py 100000 10000
$P rate.py 10000 100000 > rate.txt; $P rate.py 100000 10000 >> rate.txt

# 2. arithmetic features of the 51 exponents against the heuristic (results.txt), and the independent recheck
(cd structure && run $P analysis.py > /dev/null && run $P check.py > check.txt)

# 3. Lucas-Lehmer tests per exponent: deeper trial factoring, P-1, orderings, the ceiling
(cd touch && run $P tfdeep.py 10000 100000 10000000 && run $P tfdeep.py 100000 10000 1000000 \
  && run $P tftable.py > tftable.out && run $P tfcheck.py > tfcheck.out \
  && run $P pm1.py 10000 10000 3 10000 100,300,1000,3000,10000,30000,100000 1000 3000,10000,30000,100000,1000000 pm1_1e4.json \
  && run $P pm1.py 100000 1000000 10000 40000 100,300,1000 1000 3000 pm1_1e5a.json \
  && run $P pm1.py 100000 1000000 40000 70000 100,300,1000 1000 3000 pm1_1e5b.json \
  && run $P pm1.py 100000 1000000 70000 100000 100,300,1000 1000 3000 pm1_1e5c.json 3 \
  && run $P pm1table.py 1e4 > pm1table_1e4.out && run $P pm1table.py 1e5 > pm1table_1e5.out \
  && run $P order.py > order.out && run $P bound.py > bound.out && run $P verify.py > verify.out)

# 4. the page's engine, cut out of ../riemann/prime.html, against the list, BigInt and gmpy2; the readable engine's own test
(cd engine && run node page_check.js > page_check.json 2> page_check.txt && run $P page_check.py < page_check.json >> page_check.txt \
  && run node test_best.js > test_best.txt)

# 5. the exponents as plain numbers: successor maps, closed formulas, OEIS, the page's transformations (numbers/)
(cd numbers && run $P seqs.py 250 250 > seqs.out \
  && for s in d1 d2 d3 tt inv-d1 inv-d2 inv-d3; do run $P search.py $s > log_$s.txt; done \
  && run $P report.py && run $P formulas.py && run $P formulas_params.py && run $P formulas_oos.py > formulas_oos.txt \
  && run $P oeis.py && run node page142.js > page142.json && run $P page142_report.py && run $P verify.py > verify.txt \
  && run $P a366029.py 4500 > a366029.txt)

# 6. integer relations, bit strings, bases, growth (relations/)
(cd relations && run $P gen_nulls.py 500 \
  && run $P angle1_linrec.py > angle1.log && run $P angle2_relations.py > angle2_relations.log && run $P angle2b_tail.py \
  && run $P angle3_bits.py > angle3_bits.log && run $P angle4_bases.py > angle4_bases.log && run $P angle5_growth.py > angle5_growth.log \
  && run $P combine.py && run $P oos_select.py > oos_select.out && run $P check_real.py > check_real.out && run $P check_base23.py > check_base23.out)
