#-------------------------------------------------------# 
# Generate LAMMPS data file using random atom positions #
# 														#
# 	- Bradley Huddleston, February 27, 2018				#
#-------------------------------------------------------#

import numpy as np

# Number of atoms to create
natoms = 1000

# Size of the system cell in Angstroms
#	assuming an cubic cell, starting at the origin
system_size = 20.0

# Generate atom positions
#	Randomness for amorphous glass
positions = []
for i in range(natoms):
	positions.append(np.random.rand(3)*system_size)
	
# Write LAMMPS data file
with open('random.data','w') as fdata:
	# First line is a comment line 
	fdata.write('Random atoms - written for EnCodeVentor tutorial\n\n')
	
	#--- Header ---#
	# Specify number of atoms and atom types 
	fdata.write('{} atoms\n'.format(natoms))
	fdata.write('{} atom types\n'.format(1))
	# Specify box dimensions
	fdata.write('{} {} xlo xhi\n'.format(0.0, system_size))
	fdata.write('{} {} ylo yhi\n'.format(0.0, system_size))
	fdata.write('{} {} zlo zhi\n'.format(0.0, system_size))
	fdata.write('\n')

	# Atoms section
	fdata.write('Atoms\n\n')

	# Write each position 
	for i,pos in enumerate(positions):
		fdata.write('{} 1 {} {} {}\n'.format(i+1,*pos))
		
	# If you have bonds and angles, further sections below
	