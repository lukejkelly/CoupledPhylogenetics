#!/bin/bash

for FILE in *.supp; do
    sleep 0.25
    sbatch --export=ALL,TARGET=${FILE/.supp/} ../../src/spr_matrix.sl
done
