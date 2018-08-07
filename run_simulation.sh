#!/bin/bash

inds=($(seq 0 1 2))

for i in "${inds[@]}"; do
    python run_grow_trees_and_run_cpd.py "$i"
done


