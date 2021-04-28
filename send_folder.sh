#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Pass name of folder"
else
    rsync -ahHSv "$1" ceremade:CoupledPhylogenetics
fi
