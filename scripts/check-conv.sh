#!/bin/bash

#This script checks which directory has converged VASP relaxations
parent_directory="./"
mkdir converged
converged_directory="./converged/"

# Iterate over directories starting with "co2" in the parent directory
for dir in "$parent_directory"/co2*/; do
    if [ -f "$dir/OUTCAR" ]; then
        if grep -q "reached required accuracy" "$dir/OUTCAR"; then
            echo "Moving $dir to converged directory"
            #mv "$dir" "$converged_directory"
        fi
    fi
done

