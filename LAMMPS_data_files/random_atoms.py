#
#
#
#
#

import numpy as np

natoms = 100
cell_size = 10.0

positions = []
for i in range(natoms):
    positions.append(np.random.rand(3)*cell_size)
    

with open('main_atoms.data','w') as lammpsdata:

    lammpsdata.write('LAMMPS data file written for an EnCodeVentor tutorial\n\n')
    
    lammpsdata.write('{} atoms\n'.format(natoms))
    lammpsdata.write('1 atom types\n')
    lammpsdata.write('{} {} xlo xhi\n'.format(0.0, cell_size))
    lammpsdata.write('{} {} ylo yhi\n'.format(0.0, cell_size))
    lammpsdata.write('{} {} zlo zhi\n'.format(0.0, cell_size))
    lammpsdata.write('\n')
    
    lammpsdata.write('Atoms\n\n')
    for i, pos in enumerate(positions):
    
        lammpsdata.write('\t{} {} {} {} {}\n'.format(i+1,1,*pos))
        
    # bonds and such below
    
    

