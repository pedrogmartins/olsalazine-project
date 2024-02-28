#!/bin/bash

#This scrip goes through all relaxations in current directory and resubmits all 
#   calculations that have not converged.


for file in lammps*; do
    # Create a directory with the same name followed by -2nd
    mkdir "${file}-2nd"
    cd $file
    cp CONTCAR POTCAR run_perlmt.sh INCAR KPOINTS ../"${file}-2nd"
    mv ../"${file}-2nd"/CONTCAR ../"${file}-2nd"/POSCAR
    cd ..
    cd ./"${file}-2nd"
    pwd
    sbatch run_perlmt.sh
    cd ..
done
