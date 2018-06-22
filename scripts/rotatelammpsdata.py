#!/bin/env python
#-------------------------------------------------------# 
# Rotate atoms from a LAMMPS data file  			    #
# 														#
# 	- Bradley Huddleston, May 26, 2018					#
#-------------------------------------------------------#

"""

"""

import numpy as np
from itertools import product
import argparse
import math
import warnings


#####---------------------------------- Function definitions ----------------------------------#####
# Rotation matrices from assorted forms
def rotation_axisangle(axis, angle):
	a = math.cos(angle/2.0)
	s = math.sin(angle/2.0)
	b = axis[0]*s
	c = axis[1]*s 
	d = axis[2]*s 
	
	return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d+a*c)],
					 [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
					 [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
					 
def rotation_quaternion(q):
    
    axis = np.array(q[:3])/np.linalg.norm(q[:3])
    angle = 2*np.arctan2(np.linalg.norm(q[:3]),q[3])
    
    return rotation_axisangle(axis,angle)
    #a = q[3]
    #b = q[0]
    #c = q[1]
    #d = q[2]
    
    #return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d+a*c)],
	#				 [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
	#				 [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
					 
def rotation_eulerxzx(alpha, beta, gamma):

    s1 = np.sin(alpha)
    c1 = np.cos(alpha)
    s2 = np.sin(beta)
    c2 = np.cos(beta)
    s3 = np.sin(gamma)
    c3 = np.cos(gamma)

    return np.array([[c2,    -c3*s2,         s2*s3],
					 [c1*s2, c1*c2*c3-s1*s3, -c3*s1-c1*c2*s3],
					 [s1*s2, c1*s3+c2*c3*s1, c1*c3-c2*s1*s3]])

def write_datafile(filename, atoms, types=None, 
            xlo=0.0, ylo=0.0, zlo=0.0, xhi=1.0, yhi=1.0, zhi=1.0):
        
    # If no types, assume one type    
    if types is None:
        types = np.ones(len(atoms),dtype=int)
        
    # Write LAMMPS data file
    with open(filename,'w') as fdata:
        # First line is a comment line 
        fdata.write('Crystalline Al atoms - written for EnCodeVentor tutorial\n\n')

        #--- Header ---#
        # Specify number of atoms and atom types 
        fdata.write('{} atoms\n'.format(len(atoms)))
        fdata.write('{} atom types\n'.format(len(np.unique(types))))
        # Specify box dimensions
        fdata.write('{} {} xlo xhi\n'.format(xlo, xhi))
        fdata.write('{} {} ylo yhi\n'.format(ylo, yhi))
        fdata.write('{} {} zlo zhi\n\n'.format(zlo, zhi))

        # Atoms section
        fdata.write('Atoms\n\n')

        # Write each position 
        for i,data in enumerate(zip(atoms,types)):
            pos, t = data
            fdata.write('{} {} {} {} {}\n'.format(i+1,int(t),*pos))
	
        # If you have bonds and angles, further sections below
        
        
def read_datafile(filename):
    
    with open(filename, 'r') as f:
        # Comment line
        f.readline()
        # Empty line
        f.readline()
        
        natoms = int(f.readline().split()[0])
        ntypes = int(f.readline().split()[0])
        xlo, xhi = [float(l) for l in f.readline().split()[:2]]
        ylo, yhi = [float(l) for l in f.readline().split()[:2]]
        zlo, zhi = [float(l) for l in f.readline().split()[:2]]
        
        f.readline()
        # Atoms header
        f.readline()
        f.readline()
        
        atoms = [np.zeros(3) for i in range(natoms)]
        types = np.zeros(natoms)
        for i in range(natoms):
            line = f.readline().split()
            aid  = int(line[0])-1
            tp   = int(line[1])
            atoms[aid] = np.array([float(j) for j in line[2:5]])
            types[aid] = tp
            
    return atoms, types, xlo, ylo, zlo, xhi, yhi, zhi
            
        

#####------------------------------------------ main ------------------------------------------#####
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='rotatelammpsdata.py',
    description='''Rotate atoms in a LAMMPS data file by an 
       arbitrary rotation. The results are written to the file 
       specified by the required "out" option. The rotation 
       can be specified by (XZX) Euler angles, or an 
       axis-angle. The origin of rotation can also be 
       specified, by default it is the origin.
       
       Example: rotatelammpsdata.py --deg --axisangle 1 0 0
                30 --out rotated.data lammps.data''')
    parser.add_argument('-o','--out', metavar='FILE', type=str, required=True,
                        help='Filename for the output LAMMPS datafile')
    parser.add_argument('-d','--deg', action='store_true', default=False,
                        help='flag indicates angles are specified in degrees')
    parser.add_argument('-r','--rad', action='store_true', default=True,
                        help='flag indicates angles are specified in radians [default]')

    parser.add_argument('-O','--origin', metavar='N', nargs=3, type=float, default=[0.0, 0.0, 0.0],
                        help='Point to use as the origin of rotation')                        
                        
    # Specify rotation by axis/angle
    parser.add_argument('-a','--axisangle', metavar='N', nargs=4, type=float,
                        help='Axis/Angle representation of rotation')
    # Specify rotation by quaternion
    parser.add_argument('-q','--quat', metavar='N', nargs=4, type=float,
                        help='Quaternion representation of rotation [x y z w]')
    # Specify rotation by XZX Euler Angles
    parser.add_argument('-e','--euler', metavar='A', nargs=3, type=float,
                        help='XZX Euler angles of rotation')
    # Positional argument for filename
    parser.add_argument('infile', metavar='FILENAME', type=str, 
                        help='Input file in LAMMPS data file format')

    args = parser.parse_args()
    
    if args.deg:
        args.rad = False
        
    origin = np.array(args.origin)
    
    if args.axisangle is not None:
        axis = args.axisangle[:3]
        if args.deg:
            angle = args.axisangle[3]*np.pi/180
        else: 
            angle = args.axisangle[3]
        rot = rotation_axisangle(np.array(axis)/np.linalg.norm(axis), angle)
    elif args.quat is not None:
        rot = rotation_quaternion(args.quat)
    elif args.euler is not None:
        if args.deg:
            angles = [a*np.pi/180. for a in args.euler]
        else:
            angles = args.euler
        rot = rotation_eulerxzx(*angles)
    else:
        warnings.warn("""
Rotation must be specified either by axis/angle, quaternion, or euler angles 
from the command line. Use the -h option to see syntax help.""", SyntaxWarning)
    
    # Read in atom positions from datafile
    atoms, types, xlo, ylo, zlo, xhi, yhi, zhi = read_datafile(args.infile)
    
    # Rotate each atom around origin of rotation
    newatoms = np.array(atoms).copy()
    for i, atom in enumerate(atoms):
        newatoms[i] = np.matmul(rot, atom-origin) + origin
    
    # Find max and min in each direction for system cell
    xlo, ylo, zlo = np.min(newatoms,axis=0)
    xhi, yhi, zhi = np.max(newatoms,axis=0)    
    
    # Write lammps data file
    write_datafile(args.out, newatoms, types,
            xlo=xlo, ylo=ylo, zlo=zlo, xhi=xhi, yhi=yhi, zhi=zhi)

