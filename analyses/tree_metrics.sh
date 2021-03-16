#!/bin/bash

# run from simulations output folder
for FILE in *_x.nex; do
    sleep 0.25
    qsub -d . -v TARGET=${FILE/_x.nex/} ../../analyses/tree_metrics.pbs
done
