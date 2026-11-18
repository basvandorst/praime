#!/bin/bash
# the three-cell wheel on the published list (two halves) and alone on every odd number below 1e8 that is not prime
for part in 0 1; do bzcat ../feitsma/psps64.txt.bz2 | timeout 3000 nice ./fold -three -s 2 $part > three_list_${part}.out 2> three_list_${part}.err; echo done >> three_list_${part}.err; done
./oddcomp 3 100000000 | timeout 3000 nice ./fold -three > three_alone_1e8.out 2> three_alone_1e8.err; echo done >> three_alone_1e8.err
