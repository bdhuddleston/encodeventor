#-------------------------------------------------------# 
# Generate LAMMPS data file using random atom positions #
# 														#
# 	- Bradley Huddleston, February 27, 2018				#
#-------------------------------------------------------#

import numpy as np

# Lattice parameter for aluminum
lattice_parameter = 4.046

# Cubic FCC basis 
basis = np.array([[1.0, 0.0, 0.0],
				  [0.0, 1.0, 0.0],
				  [0.0, 0.0, 1.0]])*lattice_parameter
				  
base_atoms = np.array([[0.0, 0.0, 0.0],
					   [0.5, 0.5, 0.0],
					   [0.5, 0.0, 0.5],
					   [0.0, 0.5, 0.5]])*lattice_parameter

# Size of the system cell in lattice units
#	assuming an cubic cell starting at the origin
system_size = 20

# Generate atom positions
positions = []
for i in range(system_size):
	for j in range(system_size):
		for k in range(system_size):
			base_position = np.array([i,j,k])
			cart_position = np.inner(basis.T, base_position)
			for atom in base_atoms:
				positions.append(cart_position + atom)
	 
# Write LAMMPS data file
with open('crystalline.data','w') as fdata:
	# First line is a comment line 
	fdata.write('Crystalline Al atoms - written for EnCodeVentor tutorial\n\n')
	
	#--- Header ---#
	# Specify number of atoms and atom types 
	fdata.write('{} atoms\n'.format(len(positions)))
	fdata.write('{} atom types\n'.format(1))
	# Specify box dimensions
	fdata.write('{} {} xlo xhi\n'.format(0.0, system_size*lattice_parameter))
	fdata.write('{} {} ylo yhi\n'.format(0.0, system_size*lattice_parameter))
	fdata.write('{} {} zlo zhi\n'.format(0.0, system_size*lattice_parameter))
	fdata.write('\n')

	# Atoms section
	fdata.write('Atoms\n\n')

	# Write each position 
	for i,pos in enumerate(positions):
		fdata.write('{} 1 {} {} {}\n'.format(i+1,*pos))
		
	# If you have bonds and angles, further sections below
	