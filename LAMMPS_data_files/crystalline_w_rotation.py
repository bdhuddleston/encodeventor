#-------------------------------------------------------# 
# Generate LAMMPS data file using crystalline metal 	#
# 	atom positions 									    #
# 														#
# 	- Bradley Huddleston, May 26, 2018					#
#-------------------------------------------------------#

import numpy as np

from itertools import product

### New code - rotation matrix from axis-angle
import math

def rotation_matrix(axis, angle):
	a = math.cos(angle/2.0)
	s = math.sin(angle/2.0)
	b = axis[0]*s
	c = axis[1]*s 
	d = axis[2]*s 
	
	return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d+a*c)],
					 [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
					 [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
### end new code 

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

axis = [1, 1, 0]
angle = np.pi/6
rot = rotation_matrix(axis, angle)

basis = np.matmul(rot.T, basis)
base_atoms = np.array([np.matmul(rot,atom) for atom in base_atoms])
					   
# Size of the system cell in angstroms
#	assuming an cubic cell starting at the origin
system_size = 20.0
systembasis = np.array([[1.0, 0.0, 0.0],
						[0.0, 1.0, 0.0],
						[0.0, 0.0, 1.0]])*system_size
						
invbasis = np.linalg.inv(basis)
invsystem = np.linalg.inv(systembasis)

corners = np.array(list(product([0.0, system_size],repeat=3)))
extents = np.array([np.matmul(invbasis.T,corner) for corner in corners])
mins = [int(np.min([m,0])) for m in np.floor(np.min(extents,axis=0))]
maxs = [int(np.max([m,0]))+1 for m in np.ceil(np.max(extents,axis=0))]

# Generate atom positions
positions = []
for i in range(mins[0],maxs[0]):
	for j in range(mins[1],maxs[1]):
		for k in range(mins[2],maxs[2]):
			base_position = np.array([i,j,k])
			cart_position = np.inner(basis.T, base_position)
			for atom in base_atoms:
				atompos = cart_position + atom
				systempos = np.matmul(invsystem.T, atompos)
				if all(systempos >= 0.0) and all(systempos < 1.0):
					positions.append(atompos)
	 
# Write LAMMPS data file
with open('rotatedcrystalline.data','w') as fdata:
	# First line is a comment line 
	fdata.write('Crystalline Al atoms - written for EnCodeVentor tutorial\n\n')
	
	#--- Header ---#
	# Specify number of atoms and atom types 
	fdata.write('{} atoms\n'.format(len(positions)))
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
	