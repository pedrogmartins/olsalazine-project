#!/bin/bash


#This script identifies all simulations that have timed out and 
#  mvoes them to another directory for resubmission

# Define the root directory
root_directory="."

# Define the directory where you want to move the folders
timeout_directory="$root_directory/timeout"

# Create the timeout directory if it doesn't exist
mkdir -p "$timeout_directory"

# Iterate through all directories in the root directory
for directory in "$root_directory"/lammps*; do
    if [ -d "$directory" ]; then
        # Check if there's a *.err file containing the expression
        if grep -q "DUE TO TIME" "$directory"/*.err; then
            # Move the directory to the timeout folder
            #echo $directory
            mv "$directory" "$timeout_directory"
        fi
    fi
done
