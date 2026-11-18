#!/bin/bash
# usage: bigrun.sh PART(0|1) TAG MACHINE...   (reads every second line of the published list)
part=$1; tag=$2; shift 2
bzcat ../feitsma/psps64.txt.bz2 | awk -v p=$part 'NR%2==p' | nice ${TOURBIN:-./tour} stdin "$@" > big_${tag}_${part}.out 2> big_${tag}_${part}.err
echo done >> big_${tag}_${part}.err
