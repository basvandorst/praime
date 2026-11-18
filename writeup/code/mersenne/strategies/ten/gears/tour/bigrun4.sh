#!/bin/bash
# usage: bigrun4.sh PART(0|1) TAG RING...   (the folded-wheel program reads every second line of the published list)
part=$1; tag=$2; shift 2
bzcat ../feitsma/psps64.txt.bz2 | timeout 5400 nice ./fold -s 2 $part "$@" > big4_${tag}_${part}.out 2> big4_${tag}_${part}.err
echo done >> big4_${tag}_${part}.err
