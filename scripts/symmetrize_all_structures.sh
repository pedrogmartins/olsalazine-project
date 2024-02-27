#!/bin/bash

#This script executes symmetrize.py script in each subdirectory
#   in this folder

# Get the current directory (mother directory)
mother_directory=$(pwd)

# Loop through directories that start with "lammps"
for directory in lammps*/; do
    # Change to the lammps directory
    cd "$directory"
        
    pwd 

    # Execute the symmetrize.py script
    python3 "$mother_directory/symmetrize-P3_and_P3221.py"
    
    # Return to the mother directory
    cd "$mother_directory"
done

echo "Finished executing symmetrize-P3_and_P3221.py in lammps directories."

