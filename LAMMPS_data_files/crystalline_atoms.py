#
#
#
#
#

import numpy as np

#natoms = 100
cell_size = 4.046*3

# FCC
basis = np.array([[1.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0]])
basis_atoms = np.array([[0.0, 0.0, 0.0],
                        [0.5, 0.5, 0.0],
                        [0.5, 0.0, 0.5],
                        [0.0, 0.5, 0.5]])
lattice_constant = 4.046

positions = []
newposition = np.array([0.0, 0.0, 0.0])
simplecell = np.ones(3)*cell_size;
while((newposition[0] + 0.5*lattice_constant) < cell_size):
    while((newposition[1] + 0.5*lattice_constant) < cell_size):
        while((newposition[2] + 0.5*lattice_constant) < cell_size):
            for pos in basis_atoms:
                positions.append(newposition + pos*lattice_constant)
            newposition += basis[2]*lattice_constant
        newposition[2] = 0.0 # cell_size
        newposition += basis[1]*lattice_constant
    newposition[1] = 0.0 # cell_size
    newposition += basis[0]*lattice_constant

with open('crystalline_atoms.data','w') as lammpsdata:

    lammpsdata.write('LAMMPS data file written for an EnCodeVentor tutorial\n\n')
    
    lammpsdata.write('{} atoms\n'.format(len(positions)))
    lammpsdata.write('1 atom types\n')
    lammpsdata.write('{} {} xlo xhi\n'.format(0.0, cell_size))
    lammpsdata.write('{} {} ylo yhi\n'.format(0.0, cell_size))
    lammpsdata.write('{} {} zlo zhi\n'.format(0.0, cell_size))
    lammpsdata.write('\n')
    
    lammpsdata.write('Atoms\n\n')
    for i, pos in enumerate(positions):
    
        lammpsdata.write('\t{} {} {} {} {}\n'.format(i+1,1,*pos))
        
    # bonds and such below
    
    

