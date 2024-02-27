import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
import os
import csv

#This script goes through all folders in this directory, finds the
#   LAMMPS output file "log.lammps" and extracts the simulation box
#   parameters (converted to lattice paraemters) and energy. Both
#   are saved as .csv files for later equilibration and thermal
#   expansion analysis. This script also extract the collected
#   trajectories of all atoms in the simulation for diffusion analysis.

#Go to simulation folder
cd ../simulations

#Define needed functions

#Function to extract energy
def ext_energy(log_file):

    # Open the log.lammps file
    with open(log_file, 'r') as file:
        # Find the line with the header
        header_line = None
        found_MD = False
        #Boolean to find beginning of data output
        for line in file:
            if line.startswith("run 1000000") and found_MD == False:
                found_MD = True
            if line.startswith("Step Temp PotEng TotEng Enthalpy Fmax Fnorm") and found_MD:
                header_line = line
                break

        # If the header line is found, process the data
        if header_line is not None:
            # Find the indices of the columns we're interested in
            columns = header_line.split()
            step_index = columns.index("Step")
            temp_index = columns.index("Temp")
            poteng_index = columns.index("PotEng")

            # Initialize lists to store the data
            steps = []
            temperatures = []
            pot_energies = []

            # Read the lines and extract the data
            for line in file:
                data = line.split()
                if len(data) >= max(step_index, temp_index, poteng_index) + 1:
                    steps.append(int(data[step_index]))
                    temperatures.append(float(data[temp_index]))
                    pot_energies.append(float(data[poteng_index]))
                    if int(data[step_index]) == 1000000:
                        break
                else:
                    break

        return steps, temperatures, pot_energies

#Function to extract energy
def extract_lat_par(dat_file):

    # Initialize a dictionary to store the coordinates for each atom
    atom_coordinates = {}
    lammps_box = {}
    lammps_box[0] = {'lo_bound': [], 'hi_bound': [], 'tilt': []}
    lammps_box[1] = {'lo_bound': [], 'hi_bound': [], 'tilt': []}
    lammps_box[2] = {'lo_bound': [], 'hi_bound': [], 'tilt': []}

    steps = []

    #Equilibration step, to be editted later
    eq_step = 0

    # Open the file
    with open(dat_file, 'r') as file:
        #Booleans to orient data extraction
        atom_data_started = False
        atom_timestep_started = False
        lattice_started = False
        above_eq = False

        #Look over line
        for line in file:
            parts = line.strip().split()

            #Collect steps if above equilibration step
            if atom_timestep_started:
                step = parts[0]
                atom_timestep_started = False
                if int(step) > eq_step:
                    above_eq = True
                    steps.append(int(step))

            #Extract lattice parameters
            if lattice_started and above_eq and parts[0] != 'ITEM:':
                i, j, k = float(parts[0]), float(parts[1]), float(parts[2])
                lammps_box[counter]["lo_bound"].append(i)
                lammps_box[counter]["hi_bound"].append(j)
                lammps_box[counter]["tilt"].append(k)
                counter += 1

            #Extract coordinates of atoms
            if atom_data_started and above_eq and len(parts) >= 8:
                atom_id = int(parts[0])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])

                if atom_id not in atom_coordinates:
                    atom_coordinates[atom_id] = {'x': [], 'y': [], 'z': []}

                atom_coordinates[atom_id]['x'].append(x)
                atom_coordinates[atom_id]['y'].append(y)
                atom_coordinates[atom_id]['z'].append(z)

            #Found a new frame
            if len(parts) > 1 and parts[0] == 'ITEM:' and parts[1] == 'TIMESTEP':
                atom_data_started = False
                atom_timestep_started = True

            #Found lattice paraemter
            if len(parts) > 1 and parts[1] == 'BOX':
                lattice_started = True
                atom_timestep_started = False

            #Found atoms to extract coordinates from
            if len(parts) > 1 and parts[0] == 'ITEM:' and parts[1] == 'ATOMS':
                atom_data_started = True
                lattice_started = False
                counter = 0

    #Convert to lattice parameters
    lattice_parameters = {'a': [], 'b': [], 'c': []}
    lammps_box_2 = {'xlo': [], 'xhi': [], 'ylo': [], 'yhi': [], 'zlo': [], 'zhi': [], 'xy': [], 'xz': [], 'yz': []}


    for index in range(0, len(steps)):

        xlo_bound = lammps_box[0]["lo_bound"][index]
        xhi_bound = lammps_box[0]["hi_bound"][index]
        xy =lammps_box[0]["tilt"][index]
        ylo_bound  = lammps_box[1]["lo_bound"][index]
        yhi_bound = lammps_box[1]["hi_bound"][index]
        xz = lammps_box[1]["tilt"][index]
        zlo_bound = lammps_box[2]["lo_bound"][index]
        zhi_bound = lammps_box[2]["hi_bound"][index]
        yz = lammps_box[2]["tilt"][index]

        xlo = xlo_bound - min(0, xy, yz, xy + xz)
        xhi = xhi_bound - max(0, xy, yz, xy + xz)
        ylo = ylo_bound - min(0, yz)
        yhi = yhi_bound + max(0, yz)
        zlo = zlo_bound
        zhi = zhi_bound

        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo

        a = lx
        b = np.sqrt(ly**2 + xy**2)
        c = np.sqrt(lz**2 + xz**2 + yz**2)

        lattice_parameters["a"].append(a)
        lattice_parameters["b"].append(b)
        lattice_parameters["c"].append(c)
        lammps_box_2["xlo"].append(xlo)
        lammps_box_2["xhi"].append(xhi)
        lammps_box_2["ylo"].append(ylo)
        lammps_box_2["yhi"].append(yhi)
        lammps_box_2["zlo"].append(zlo)
        lammps_box_2["zhi"].append(zhi)
        lammps_box_2["xy"].append(xy)
        lammps_box_2["xz"].append(xz)
        lammps_box_2["yz"].append(yz)

    return lattice_parameters

#Extract data from co2 bound simulations

# Define the base directory where your folders are located
base_dir = "./"

# Define the folder names
folders = ["co2-0.1K", "co2-1K", "co2-10K", "co2-25K", "co2-50K", "co2-100K", "co2-150K", "co2-200K", "co2-250K"]

#Extract energies and lattice parameters
energy_arrays = []
lat_par_dictionaries = []

# Iterate through each folder
for folder in folders:
    folder_path = os.path.join(base_dir, folder)
    # Check if the folder exists
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        print(f"Files ending with '.dat' in folder {folder}:")
        # List files in the folder
        for file in os.listdir(folder_path):
            if file.endswith("log.lammps"):
            #Extract energies
                a, b, c = ext_energy(os.path.join(folder_path, file))
                energy_arrays.append(c)
            if file.endswith(".dat"):
                #Extract lattice parameters
                lat_par_dictionaries.append(extract_lat_par(os.path.join(folder_path, file)))

#Save data as .csv for later analysis 
f = open('co2_md-data_lat-par.csv','wb')
w = csv.DictWriter(f,lat_par_dictionaries)
f.close()

f = open('co2_md-data_energy.csv','wb')
w = csv.DictWriter(f,energy_arrays)
f.close()
