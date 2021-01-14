#!/bin/bash

for FILE in *_x.nex; do
    sleep 0.25
    qsub -d . -v TARGET=${FILE/_x.nex/} ../../src/tree_metrics.pbs
done
