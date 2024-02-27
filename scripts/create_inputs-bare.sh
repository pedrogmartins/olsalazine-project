#!/bin/bas

for i in $(ls best-30-structures/ | grep bare); do
   mkdir ${i:0:-5}
   cp best-30-structures/$i ${i:0:-5}/POSCAR
   cp base-files/* ${i:0:-5}
   cd ${i:0:-5}
   sbatch run_perlmt.sh
   cd ..
done
