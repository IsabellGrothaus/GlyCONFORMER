#/bin/bash

HILLS=$1
STRIDE=$2
kt=$3

plumed sum_hills --hills $HILLS --stride $STRIDE --kt $kt

