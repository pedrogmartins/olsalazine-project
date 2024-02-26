#!/bin/bash
# script to create submit final single point calculation at the end of VASP relaxation

mkdir single_point
cp CONTCAR POTCAR run_perlmt.sh INCAR KPOINTS ./single_point
mv ./single_point/CONTCAR single_point/POSCAR
sed -i 's/IBRION = 2/IBRION = -1/g' ./single_point/INCAR
cd single_point
sed -i 's/NSW = 2000/NSW = 1/g' INCAR
sed -i 's/-q regular/-q debug/g' run_perlmt.sh
sed -i 's/-t 12:00:00/-t 00:30:00/g' run_perlmt.sh
sed -i 's/cu-mfu-4l/olz-proj/g' run_perlmt.sh
sed -i 's/pedrogm@berkeley.edu/pedrogm.jobs@gmail.com/g' run_perlmt.sh
#sbatch run_perlmt.sh



