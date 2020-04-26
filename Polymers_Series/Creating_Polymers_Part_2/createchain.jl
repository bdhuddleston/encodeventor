#
# Create a polymer chain that winds through space
#

using GeometryTypes, LinearAlgebra
using Distributions, Quaternions
const Pt3D = Point{3,Float32}

# c-c-c angle
cca = (180-109.5)*pi/180 /2

R = [cos(cca) sin(cca) 0;
    -sin(cca) cos(cca) 0;
     0        0        1]'

# mirror matrix - for hydrogens
R90 = [0  1  0;
      -1  0  0;
       0  0  1]'

# hydrogen bond rotation out of plane
Rh = [1 0 0;
      0 cos(pi/3) sin(pi/3);
      0 -sin(pi/3) cos(pi/3)]


# c-c bond length
ccl0 = 2*0.77f0
# c-h bond length
chl0 = 0.77f0 + 0.33f0

atoms = Pt3D[]
types = Int[]
bonds = NTuple{3,Int}[]     # two bond types, C-C, C-H
angles = NTuple{4,Int}[]    # two angle types, C-C-C, H-C-H
dihedrals = NTuple{4,Int}[]

# initial chain vector direction
chainvec = Pt3D(1,0,0)

# define a function that gives a random perturbation of that vector direction
function perturbate(vector, anglerange=15.0)
    # define a quaternion with a random axis and a small angle
    angle = rand(Distributions.Normal(0,anglerange/3*pi/180))

    theta = acos(2*rand() - 1)
    phi   = 2*pi*rand()
    qaxis  = [cos(theta)*cos(phi),cos(theta)*sin(phi),sin(theta)]
    q = Quaternion(cos(angle/2),sin(angle/2) .* qaxis)
    vnew = q * Quaternion([vector...]) * q'
    return axis(vnew)
end

function perturbate_uniform(vector, anglerange=15.0)
    # define a quaternion with a random axis and a small angle
    angle = (rand()*anglerange-anglerange/2)*pi/180

    theta = acos(2*rand() - 1)
    phi   = 2*pi*rand()
    qaxis  = [cos(theta)*cos(phi),cos(theta)*sin(phi),sin(theta)]
    q = Quaternion(cos(angle/2),sin(angle/2) .* qaxis)
    vnew = q * Quaternion([vector...]) * q'
    return axis(vnew)
end

function createCH2!(atoms, types, angles, dihedrals, point, vector)
    n = length(atoms)
    # create the carbon
    push!(atoms,point)
    push!(types,1)
    if n >=3
        push!(bonds,(1,n-2,n+1))
    end
    if n >= 6
        push!(angles,(1,n-5,n-2,n+1))
    end
    if n >= 9
        push!(dihedrals,(n-8,n-5,n-2,n+1))
    end

    # create the hydrogens
    push!(atoms,point + chl0 * Rh  * R90 * vector)
    push!(types,2)
    push!(bonds,(2,n+1,n+2))

    push!(atoms,point + chl0 * Rh' * R90 * vector)
    push!(types,2)
    push!(bonds,(2,n+1,n+3))
    push!(angles,(2,n+2,n+1,n+3))
    return nothing
end

while length(atoms) < 300
    global chainvec, R, R90
    R = R'
    R90 = R90'
    newcarbon = Pt3D(0,0,0)
    if length(atoms) > 2
        newcarbon = atoms[end-2] + ccl0 * R * chainvec
    end
    createCH2!(atoms, types, angles, dihedrals, newcarbon, chainvec)
    chainvec = perturbate_uniform(chainvec,90.0)
end
#push!(bonds,(1,1,length(atoms)-2))
sysx = ccl0*cos(cca)*2 * 50

open("polymerchain.data","w") do f
    print(f,"# polyethylene chain - 10 carbons \n\n")
    print(f,"\t$(length(atoms)) atoms\n") #    30 atoms
    print(f,"\t$(length(bonds)) bonds\n")
    print(f,"\t$(length(angles)) angles\n")
    print(f,"\t$(length(dihedrals)) dihedrals\n\n")
    print(f,"\t2 atom types\n")
    print(f,"\t2 bond types\n")
    print(f,"\t2 angle types\n")
    print(f,"\t1 dihedral types\n\n")
    print(f,"0.0 $sysx xlo xhi\n")
    print(f,"-30.0 30.0 ylo yhi\n")
    print(f,"-30.0 30.0 zlo zhi\n\n")
    print(f,"Masses\n\n")
    print(f,"\t1 12.011\n")
    print(f,"\t2 1.008\n\n")
    print(f,"Atoms # molecular\n\n")
    i = 0
    for (a,t) in zip(atoms,types)
        i += 1
        print(f,"\t$i 1 $t $(a[1]) $(a[2]) $(a[3])\n")
    end
    print(f,"\nBonds\n\n")
    for i in 1:length(bonds)
        t,x,y = bonds[i]
        print(f,"\t$i $t $x $y\n")
    end
    print(f,"\nAngles\n\n")
    for i in 1:length(angles)
        t,x,y,z = angles[i]
        print(f,"\t$i $t $x $y $z\n")
    end
    print(f,"\nDihedrals\n\n")
    for i in 1:length(dihedrals)
        w,x,y,z = dihedrals[i]
        print(f,"\t$i 1 $w $x $y $z\n")
    end
end
