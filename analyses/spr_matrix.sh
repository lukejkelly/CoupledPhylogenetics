#!/bin/bash

for FILE in *.supp; do
    sleep 0.25
    qsub -d . -v TARGET=${FILE/.supp/} ../../src/spr_matrix.pbs
done
