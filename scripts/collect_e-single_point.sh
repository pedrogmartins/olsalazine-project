#!/bin/bash


#This script walks through folder with converged calculations and extracts the 
# final VASP single-point calculation energy and the folder name, and stores 
# in .csv file for further processing. 

for i in $(find -type d -name "co2*"); do
   cd $i
   pwd
   cd single_point
   python3 ../../v6grad2.py OUTCAR |  grep 'Energy:' | tail -n 1 | awk '/Energy:/ {print $3}'
   cd ..
   cd ..
done > co2-worse-starting.csv
