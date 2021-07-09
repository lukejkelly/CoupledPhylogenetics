#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Pass name of folder"
else
    rsync -ahHSv ceremade:CoupledPhylogenetics/"$1"/figs "$1"
fi
