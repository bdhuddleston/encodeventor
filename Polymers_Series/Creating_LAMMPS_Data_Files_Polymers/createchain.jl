# Write polymer chain data file
using GeometryTypes, LinearAlgebra
const Pt3D = Point{3,Float32}

# start with a three carbon link

# c-c-c angle
cca = (180-109.5)*pi/180 /2

R = [cos(cca) sin(cca) 0;
    -sin(cca) cos(cca) 0;
     0        0        1]'

R90 = [0  1  0;
      -1  0  0;
       0  0  1]'

Rh = [1 0 0;
      0 cos(pi/3) sin(pi/3);
      0 -sin(pi/3) cos(pi/3)]


# c-c bond length
ccl0 = 1.33f0
chl0 = 0.99f0

atoms = Pt3D[]
types = Int[]
bonds = NTuple{3,Int}[]
angles = NTuple{3,Int}[]
dihedrals = NTuple{4,Int}[]

chainvec = Pt3D(1,0,0)

push!(atoms,Pt3D(0,0,0))
push!(types,1)

push!(atoms,Pt3D(0,chl0*cos(pi/3),chl0*sin(pi/3)))
push!(types,2)
push!(bonds,(2,1,2))

push!(atoms,Pt3D(0,chl0*cos(pi/3),-chl0*sin(pi/3)))
push!(types,2)
push!(bonds,(2,1,3))
push!(angles,(2,1,3))

i = 3
while length(atoms) < 30
    global i, chainvec, R, R90
    i += 1
    R = R'
    R90 = R90'
    push!(atoms,atoms[end-2] + ccl0 * R * chainvec)
    push!(types,1)
    push!(bonds,(1,i-3,i))
    if i > 6
        push!(angles,(i-6,i-3,i))
    end
    if i > 9
        push!(dihedrals,(i-9,i-6,i-3,i))
    end
    i += 1
    push!(atoms,atoms[end]   + chl0 * Rh  * R90 * chainvec)
    push!(types,2)
    push!(bonds,(2,i-1,i))

    i += 1
    push!(atoms,atoms[end-1] + chl0 * Rh' * R90 * chainvec)
    push!(types,2)
    push!(bonds,(2,i-2,i))
    push!(angles,(i-1,i-2,i))
end
push!(bonds,(1,1,length(atoms)-2))
sysx = ccl0*cos(cca)*2 * 5

open("polymerchain.data","w") do f
    print(f,"# polyethylene chain - 10 carbons \n\n")
    print(f,"\t$(length(atoms)) atoms\n")
    print(f,"\t$(length(bonds)) bonds\n")
    print(f,"\t$(length(angles)) angles\n")
    print(f,"\t$(length(dihedrals)) dihedrals\n\n")
    print(f,"\t2 atom types\n")
    print(f,"\t2 bond types\n")
    print(f,"\t1 angle types\n")
    print(f,"\t1 dihedral types\n\n")
    print(f,"0.0 $sysx xlo xhi\n")
    print(f,"-3.0 3.0 ylo yhi\n")
    print(f,"-3.0 3.0 zlo zhi\n\n")
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
        x,y,z = angles[i]
        print(f,"\t$i 1 $x $y $z\n")
    end
    print(f,"\nDihedrals\n\n")
    for i in 1:length(dihedrals)
        w,x,y,z = dihedrals[i]
        print(f,"\t$i 1 $w $x $y $z\n")
    end
end
