import os
import numpy as np
import shutil
import ase.io

%cd olsalz-project/phase-IV-simulated-annealing/2023_08_28-NPP-relax/bare-fully-relaxed-water-6-v1

#This script walks through a directory, finds the potential
#  energy of each structure from a single-point LAMMPS
#  computation powered by a neural network potential. It ultimately
#  prints out the folder corresponding to the 30 lowest energy
#  structures and check their number of atoms to certify structure
#  integrity.    

base_directory = './'

# List to store visited folder names
visited_folders = []
energies = []

# Loop over subdirectories
for root, dirs, files in os.walk(base_directory):
    if 'log.lammps' in files:
        folder_name = os.path.basename(root)
        visited_folders.append(folder_name)
        
        #Convert output file to vasp

        log_file_path = os.path.join(root, 'log.lammps')

        with open(log_file_path, 'r') as file:
            lines = file.readlines()

        # Define a flag to track occurrences
        
        occurrence_count = 0

        # Iterate through the lines
        for i, line in enumerate(lines):
            if "Step Temp E_pair E_mol TotEng Press" in line:
                occurrence_count += 1
                #print(folder_name)
                #print(occurrence_count)
                if occurrence_count == 3:
                    desired_line_index = i + 1
                    if desired_line_index < len(lines):
                        desired_line = lines[desired_line_index]
                        #print(desired_line)
                        # Split the line and extract the value
                        values = desired_line.split()
                        if len(values) >= 6:
                            desired_value = values[2]  # Assuming
                            #print(desired_value)

        energies.append(float(desired_value))


energies = np.asarray(energies)
visited_folders = np.asarray(visited_folders)

#Sort structures according to energy 
sorted_indices = np.argsort(energies)
sorted_energies = energies[sorted_indices]
sorted_visited_folders = visited_folders[sorted_indices]

#Extract 30 lowest energy structures
destination_directory = "best_30_structures"

for root, dirs, files in os.walk(base_directory):
    for directory in dirs:
        if directory in sorted_visited_folders[0:30]:
            dir_path = os.path.join(root, directory)
            file_to_copy = os.path.join(dir_path, "annealed_run_0.vasp")
            print(file_to_copy)
            if os.path.exists(file_to_copy):
                shutil.copy(file_to_copy, os.path.join(destination_directory, f"{directory}.vasp"))

current_directory = os.getcwd()

%cd {destination_directory}

#Double check number of atoms in each calculation to make subdirectories
#  they are sound and no atoms were lost. 

dir_path = os.path.join(current_directory, destination_directory)
   
for root, dirs, files in os.walk(dir_path):
    for file in files:
        if "bare" in file:
            print(file)
            atoms = ase.io.read(file, format = 'vasp')
            print(len(atoms))
