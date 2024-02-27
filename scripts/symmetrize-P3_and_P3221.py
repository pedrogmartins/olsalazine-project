import ase.io
from ase import Atoms
from ase.io import vasp
import ase.build
import numpy as np
import pymatgen.core as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Lattice, Structure, Element, SymmOp
from pymatgen.symmetry.groups import sg_symbol_from_int_number
import pymatgen
import sys
from ase.io.lammpsdata import write_lammps_data
import re
import os


#This script extracts the final structure in a LAMMPS simulated annealing
#  run (molecular dynamics at 400K followed by slow cooling). Then, the amines
#  (and potentially the bound water molecule) are extracted to create
#  symmetrized structures following P3 and P3321 groups symmetries for
#  VASP relaxation. This is part of a conformational search initiative.

#Define needed functions

#Confine atoms outside the unit cell with periodic boundary conditions
def to_unit_cell(fractional_coords):
    coords = []
    for i in range(len(fractional_coords)):
        a = fractional_coords[i]
        if a < 0:
            a = a + np.floor(1 - a)
        elif a > 0:
            a = a - np.ceil(a - 1)
        coords.append(a)
    return np.array(coords)

#Place atoms according to symmetr
def add_atoms(atoms, s, f):
    #print(atoms)
    ns = s.copy()
    for a in atoms:
        ns.append(a[0], a[1])
    nf = SpacegroupAnalyzer(ns, f._symprec)
    #print(nf)
    #print(len(ns))
    #print(len(s))
    ns.to(fmt='poscar', filename='temp-sym-test.vasp')
    s.to(fmt='poscar', filename='temp-test.vasp')
    #print('Previous space group: {}'.format(f.get_space_group_symbol()))
    #print('New space group: {}'.format(nf.get_space_group_symbol()))
    return ns, nf

#Sort atoms by chemical identity to match LAMMSP imput file
def sort_atoms_by_number(atoms):
    return sorted(atoms, key=lambda atom: atom.number)

current_directory = os.getcwd()
folder = os.path.basename(current_directory)
directory_name = os.path.basename(current_directory).split('-', 1)[1]

#Get name of the trajectory file from current directory and get
#   index of final structure

# Remove suffix from LAMMPS output trajectory
temp = re.search(r'-(\d+)_seed', directory_name).group(1)
if directory_name[-2:].isdigit():
    final_trajec = directory_name[:-12] + '-t' + temp + '-final_min.dat'
else:
    final_trajec = directory_name[:-11] + '-t' + temp + '-final_min.dat'

with open(final_trajec, 'r') as file:
    lines = file.readlines()

# Find the last occurrence of "ITEM: TIMESTEP"
last_timestep_index = -1
for i in range(len(lines)):
    if "ITEM: TIMESTEP" in lines[i]:
        last_timestep_index = i

if last_timestep_index != -1 and last_timestep_index + 1 < len(lines):
    index = int(lines[last_timestep_index + 1].strip())

last_structure = final_trajec[:-4] + '-' + str(index)
output_name = 'annealed_run_' + str(index) + '.vasp'

#Get the last structure from the dump file
dump_filename = final_trajec
atoms = ase.io.read(dump_filename, format='lammps-dump-text')

ase.io.lammpsdata.write_lammps_data('annealed_run_' + str(index) + '.dat', atoms)
ase.io.vasp.write_vasp('sorted-' + output_name, atoms, direct = True)


#Get all indices of the amines depending on the version of the original
#  starting structure.
v = re.search(r'v\d+', directory_name).group(0)
if v == 'v1':
    amine_1 = np.asarray([
    240, 198, 246, 216, 54, 252, 114, 108, 192, 96, 24, 162, 102, 180,
    60, 66, 186, 78, 72, 174, 90, 84, 168, 42, 48, 210, 36, 30
    ])
elif v == 'v2':
    amine_1 = np.asarray([246, 198, 252, 216, 66, 36, 192, 108, 240, 24, 30,
    162, 114, 180, 78, 72, 186, 90, 84, 174, 102, 96, 168, 54, 60, 210, 48, 42
    ])
elif v == 'v3':
    amine_1 = np.asarray([246, 198, 252, 216, 66, 36, 192, 108, 240, 24, 30,
    162, 114, 180, 78, 72, 186, 90, 84, 174, 102, 96, 168, 54, 60, 210, 48, 42
    ])
elif v == 'v4':
    amine_1 = np.asarray([240, 198, 246, 216, 54, 252, 114, 108, 192, 96, 24,
    162, 102, 180, 60, 66, 186, 78, 72, 174, 90, 84, 168, 42, 48, 210, 36, 30
    ])

amine_2 = amine_1 + 1
amine_3 = amine_1 + 2
amine_4 = amine_1 + 3
amine_5 = amine_1 + 4
amine_6 = amine_1 + 5
all_amines_single_list = np.hstack((
                                    amine_1,
                                    amine_2,
                                    amine_3,
                                    amine_4,
                                    amine_5,
                                    amine_6
                                    ))
list_of_amines = [amine_1, amine_2, amine_3, amine_4, amine_5, amine_6]

#Select the framework only
atoms = ase.io.read('sorted-' + output_name, format='vasp')
#framework_atoms = Atoms(cell = atoms.cell, pbc = atoms.pbc)

framework_atoms = Atoms(cell = atoms.cell, pbc = atoms.pbc)

for i, atom in enumerate(atoms):
    if i not in all_amines_single_list:
        framework_atoms.append(atom)

ase.io.write('framework.vasp', framework_atoms, format='vasp')

#Start arrays to store the symmetrized structures
struct_list_P3321 = []
struct_list_P3 = []

frameworks_dir = 'phase-IV-simulated-annealing/2023_09_12-all_4_structures/frameworks'


#Create structures with only one amine (P3221 symmetry)
for i in range(0, 6):
    framework = ase.io.read(frameworks_dir + '/co2-framework.vasp', format='vasp')
    atoms_to_add = Atoms([atoms[j] for j in list_of_amines[i]])
    for atom in atoms_to_add:
        framework.append(atom)
    struct_list_P3321.append(framework)

    # Sort atoms by chemical identity
    sorted_atoms = sort_atoms_by_number(framework)
    sorted_structure = Atoms(sorted_atoms)
    sorted_structure.set_cell(framework.get_cell())
    sorted_structure.set_pbc(framework.get_pbc())

    #print("P3321 symmetry")
    #print(len(framework))

    ase.io.write('amine_' + str(i) + '.vasp', sorted_structure, format='vasp')


#Create structures with two amines (P3 symmetry)
for i in range(0, 3):
    #framework = ase.io.read("framework.vasp", format='vasp')
    framework = ase.io.read(frameworks_dir + '/co2-framework.vasp', format='vasp')
    atoms_to_add = Atoms([atoms[j] for j in list_of_amines[i]])
    for atom in atoms_to_add:
        framework.append(atom)
    atoms_to_add_2 = Atoms([atoms[k] for k in list_of_amines[i+3]])
    for atom in atoms_to_add_2:
        framework.append(atom)

    struct_list_P3.append(framework)

    # Sort atoms by chemical identity
    sorted_atoms = sort_atoms_by_number(framework)
    sorted_structure = Atoms(sorted_atoms)
    sorted_structure.set_cell(framework.get_cell())
    sorted_structure.set_pbc(framework.get_pbc())

    #print("P3 symmetry")
    #print(len(framework))

    ase.io.write('amine_' + str(i) + str(i+3) + '.vasp', sorted_structure, format='vasp')

#Creating P3221 final symmetrized structures
# counter = 0
# for structure in struct_list_P3321:

#     print(counter)

#     #Get actually symmetric framework
#     sym_framework = mg.Structure.from_file(frameworks_dir + '/co2-framework.vasp')

#     # Load structures in pymatgen
#     one_amine = mg.Structure.from_file('amine_' + str(counter) + '.vasp')
#     #framework = mg.Structure.from_file("framework.vasp")
#     framework = mg.Structure.from_file(frameworks_dir + '/co2-framework.vasp')

#     #Find difference between structures
#     diff = [a for a in one_amine if a not in framework]
#     atoms = [(a.as_dict()['species'][0]['element'],a.coords) for a in diff]

#     #Get space group of actually symmetric structure
#     symprec = 0.3 #Had to increase to be able to find spacegroup
#     space_group_sym = SpacegroupAnalyzer(sym_framework, symprec)
#     print('Found spacegroup {}'.format(space_group_sym.get_space_group_symbol()))

#     #Get space group or original structure
#     space_group = SpacegroupAnalyzer(framework, symprec)
#     print('Found spacegroup {}'.format(space_group.get_space_group_symbol()))

#     sym_ops = space_group_sym.get_symmetry_operations()
#     #Loop over atoms and find equivalent sites as well as add atoms
#     for atom in atoms:
#         lattice = framework.lattice
#         coords = lattice.get_fractional_coords(atom[1])
#         equiv_sites = []
#         equiv_sites = [op.operate(coords) for op in sym_ops]
#         equiv_sites = [to_unit_cell(c) for c in equiv_sites]
#         framework, space_group_sym = add_atoms([(atom[0],e) for e in equiv_sites], framework, space_group_sym)

#     one_amine = framework.get_sorted_structure()

#     file = os.path.basename(current_directory)

#     #print(space_group.get_space_group_symbol())
#     output_folder = 'phase-IV-simulated-annealing/2023_09_12-symmetrization_II-straight-framework/P3221/co2/'
#     one_amine.to(fmt = 'poscar', filename = output_folder + file + '-all_amines_sym_' + str(counter) + '.vasp')

#     counter += 1

#Creating P3 final symmetrized structures
counter = 0
for structure in struct_list_P3:

    print(counter)

    #Get actually symmetric framework
    sym_framework = mg.Structure.from_file(frameworks_dir + '/co2-framework.vasp')

    # Symmetrize structure, load them from vasp using pymatgen

    # Load structures in pymatgen
    two_amines = mg.Structure.from_file('amine_' + str(counter) + str(counter+3) + '.vasp')
    #framework = mg.Structure.from_file("framework.vasp")
    framework = mg.Structure.from_file(frameworks_dir + '/co2-framework.vasp')

    #Find difference between structures
    diff = [a for a in two_amines if a not in framework]
    atoms = [(a.as_dict()['species'][0]['element'],a.coords) for a in diff]

    example = mg.Structure.from_file("/phase-IV-simulated-annealing/2023_09_12-symmetrization_II-straight-framework/Ba3Si6N4O9.vasp")
    symprec = 0.3 #Had to increase to be able to find spacegroup
    space_group_sym = SpacegroupAnalyzer(example, symprec)
    print('Found spacegroup {}'.format(space_group_sym.get_space_group_symbol()))

    sym_ops = space_group_sym.get_symmetry_operations()
    #Loop over atoms and find equivalent sites as well as add atoms
    for atom in atoms:
        lattice = framework.lattice
        coords = lattice.get_fractional_coords(atom[1])
        equiv_sites = []
        equiv_sites = [op.operate(coords) for op in sym_ops]
        equiv_sites = [to_unit_cell(c) for c in equiv_sites]
        framework, space_group_sym = add_atoms([(atom[0],e) for e in equiv_sites], framework, space_group_sym)

    two_amines = framework.get_sorted_structure()

    file = os.path.basename(current_directory)

    #print(space_group.get_space_group_symbol())
    output_folder = 'phase-IV-simulated-annealing/2023_09_12-symmetrization_II-straight-framework/P3/co2+water/'
    two_amines.to(fmt = 'poscar', filename = output_folder + file + '-all_amines_sym_' + str(counter) + str(counter+3) + '.vasp')

    counter += 1

