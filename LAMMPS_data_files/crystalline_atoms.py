#------------------------------------------------------------# 
# Generate LAMMPS data file using crystalline atom positions #
# 													    	 #
# 	- Bradley Huddleston, March 21, 2018				     #
#------------------------------------------------------------#

import numpy as np
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

# Decide on a lattice constant - what material?
lattice_constant = 4.046

# Pick a cell size
cell_size = 25.

cell_basis = np.array([[cell_size, 0.0, 0.0],
					   [0.0, cell_size, 0.0],
					   [0.0, 0.0, cell_size]])

# Pick a set of basis vectors (using primitive basis for this one)
#origbasis = np.array([[1.0, 0.0, 0.0],
#                  [0.0, 1.0, 0.0],
#                  [0.0, 0.0, 1.0]])
origbasis = np.array([[0.5, 0.5, 0.0],
                  [0.0, 0.5, 0.5],
                  [0.5, 0.0, 0.5]])
				  
# Create a rotation matrix to rotate [111] direction to the x-axis
angle = np.arccos(np.dot([1,0,0],[1,1,1])/np.sqrt(3))
axis = -np.cross([1,0,0],[1,1,1])
axis = axis / np.linalg.norm(axis)
rot1 = rotation_matrix(axis, angle)

# Rotate around the x-axis to align rows of atoms 
axis = [1, 0, 0]
angle = np.arccos(np.dot(np.inner(rot1.T, [np.sqrt(2)/2., np.sqrt(2.)/2., 0.0]),
					[0,1,0]))
rot2 = rotation_matrix(axis, angle)

# Combine the rotation matrices
rot = np.inner(rot2.T, rot1).T

# Rotate basis 
basis = np.inner(rot.T,origbasis).T

# Scale basis by lattice_constant
basis = basis*lattice_constant

# Atoms in the basis 
origatom = np.array([[0,0,0]]) #,
					 #[0.5, 0.5, 0.0],
					 #[0.5, 0.0, 0.5],
					 #[0.0, 0.5, 0.5]])*lattice_constant
rel_atom = np.array([np.inner(rot.T,atom) for atom in origatom])
start_position = np.array([0,0,0])

# Get furthest extent of system along each lattice vector
# Define vectors to the corners
system_corners = np.array([[0.0, 0.0, cell_size],
						  [0.0, cell_size, 0.0],
						  [cell_size, 0.0, 0.0],
						  [0.0, cell_size, cell_size],
						  [cell_size, 0.0, cell_size],
						  [cell_size, cell_size, 0.0],
						  [cell_size, cell_size, cell_size]])
# Pre-calculate inverted bases 
invbasis = np.linalg.inv(basis)
invcell_basis = np.linalg.inv(cell_basis)
# Get the extents - mins and maxs
extents = np.array([np.inner(invbasis.T, corner) for corner in system_corners])
mins = [int(np.min([m,0])) for m in np.floor(np.min(extents,axis=0))]
maxs = [int(np.max([m,0])) for m in np.ceil(np.max(extents,axis=0))]

# Iterate through ranges on basis creating positions
positions = []
for i in  range(mins[0], maxs[0]+1):
	for j in range(mins[1], maxs[1]+1):
		for k in range(mins[2], maxs[2]+1):
			bravais_point = np.array([i,j,k])
			for atom in rel_atom:
				cart_point = np.inner(basis.T,bravais_point) + atom + start_position
				cell_point = np.round(np.inner(invcell_basis.T, cart_point),6)
				if all(cell_point >= 0.0) and all(cell_point < 1.0):
					positions.append(cart_point)
				
				
# Write LAMMPS data file
with open('crystalline.data','w') as fdata:
	# First line is a comment line 
	fdata.write('Crystalline atoms - written for EnCodeVentor tutorial\n\n')
	
	#--- Header ---#
	# Specify number of atoms and atom types 
	fdata.write('{} atoms\n'.format(len(positions)))
	fdata.write('{} atom types\n'.format(1))
	# Specify box dimensions
	fdata.write('{} {} xlo xhi\n'.format(0.0, cell_size))
	fdata.write('{} {} ylo yhi\n'.format(0.0, cell_size))
	fdata.write('{} {} zlo zhi\n'.format(0.0, cell_size))
	fdata.write('\n')

	# Atoms section
	fdata.write('Atoms\n\n')

	# Write each position 
	for i,pos in enumerate(positions):
		fdata.write('{} 1 {} {} {}\n'.format(i+1,*pos))
		
	# If you have bonds and angles, further sections below

    

