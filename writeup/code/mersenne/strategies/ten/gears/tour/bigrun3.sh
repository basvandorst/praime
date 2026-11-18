#!/bin/bash
# usage: bigrun3.sh PART(0|1) TAG RING...   (ring3 reads every second line of the published list)
part=$1; tag=$2; shift 2
bzcat ../feitsma/psps64.txt.bz2 | timeout 5400 nice ./ring3 -s 2 $part "$@" > big3_${tag}_${part}.out 2> big3_${tag}_${part}.err
echo done >> big3_${tag}_${part}.err
