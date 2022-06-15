import sys
import os
from ase import io, Atoms
from ase.data import atomic_numbers, atomic_masses
from numpy import *
from math import sqrt

basepath=os.getcwd()

#pseudodir='/global/common/cori/software/vasp/pseudopotentials/PBE/potpaw_PBE.54/'

str_name=sys.argv[1]

def read_POSCAR(str_name):
    infile=os.path.join('%s.vasp')%(str_name)
    infile_exists=os.path.exists(infile)
    if not infile_exists:
        print("POSCAR file does not exists")
    rf = open(infile,'rt')
    comment_line=rf.readline()
    lattice_scale=float(rf.readline().split()[0])
    lattice_vectors=[]
    for line in range(3):
        lat_abc=rf.readline().split()
        floatvec=float(lat_abc[0]),float(lat_abc[1]),float(lat_abc[2])
        lattice_vectors.append(floatvec)

    basis_vectors=array(lattice_vectors)*lattice_scale
    cellX=basis_vectors[0]
    cellY=basis_vectors[1]
    cellZ=basis_vectors[2]
    cell=(cellX, cellY, cellZ)

    atom_symbols=[]
    numofatoms=rf.readline().split()
    atomtypes=numofatoms
    numofatoms=rf.readline().split()
    numsys = len(numofatoms)
    for i, num in enumerate(numofatoms):
        numofatoms[i] = int(num)
        [atom_symbols.append(atomtypes[i]) for na in range(numofatoms[i])]
    coord_sys=rf.readline()
    cartesian = coord_sys[0].lower() == 'c' or coord_sys[0].upper() == 'C'
    tot_natoms=sum(numofatoms)
    atoms_pos=empty((tot_natoms,3))
    for atom in range(tot_natoms):
        ap = rf.readline().split()
        atoms_pos[atom]=(float(ap[0]),float(ap[1]),float(ap[2]))
    if cartesian:
        atoms_pos *= lattice_scale
    atoms = Atoms(symbols=atom_symbols, cell=basis_vectors, pbc=True)
    if cartesian:
        atoms.set_positions(atoms_pos)
    else:
        atoms.set_scaled_positions(atoms_pos)
    #print(atoms,lattice_scale,cellX,cellY,cellZ)
    return atoms,lattice_scale,cellX,cellY,cellZ
    

def create_POSCAR(atoms,lattice_scale,cellX,cellY,cellZ,str_name,cartesian):
    if cartesian == 'True' :
        coord=atoms.get_positions()
    else:
        coord=atoms.get_scaled_positions()

    old_cell=(cellX,cellY,cellZ)

    print(old_cell)

    cella_hcp=sqrt(sum(cellX**2))
    cellb_hcp=sqrt(sum(cellY**2))
    cellc_hcp=sqrt(sum(cellZ**2))
    
    volume_hcp= 3*sqrt(2)*cella_hcp**3


    transform_matrix= array([[(1./2.),0.,0.],[0.,(sqrt(3)/2.),0.],[0,0,1]]) 

# new cell generation

    cella_ortho,cellb_ortho,cellc_ortho=dot(transform_matrix.transpose(),(cella_hcp,cellb_hcp,cellc_hcp))

    new_cell=[cella_ortho,cellb_ortho,cellc_ortho,90,90,90]

    volume_ortho= cella_ortho*cellb_ortho*cellc_ortho

    frac_vol=(volume_ortho/volume_hcp)

    print(new_cell)

    P_matrix=linalg.inv(transform_matrix.transpose())
 
    print(P_matrix)

    old_origin=[0.,0.,0.]
    
# new origin

    new_origin=dot(P_matrix,old_origin)

    print(new_origin)
	
    atoms.set_cell(new_cell)

# new coordinate

    atoms.wrap()   
  
    ind_sort = argsort(atoms.get_chemical_symbols())
    symbols = array(atoms.get_chemical_symbols())[ind_sort]
    new_coord = coord[ind_sort]
    print(new_coord)

    outfile='%s_ortho.vasp'%(str_name)

#writing a description of the output file
    new_comment="hcp to orthorhombic"
    
    of = open(outfile,'w')
    of.write(new_comment + '\n')

#writing the basis vectors for the new cell

    long_format=True

    of.write('%19.16f\n' % lattice_scale)
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'
    for vec in atoms.get_cell():
        of.write(' ')
        for el in vec:
            of.write(latt_form % el)
        of.write('\n')

#writing the atom symbols
    sc = []
    psym = symbols[0]
    count = 0
    for sym in symbols:
        if sym != psym:
            sc.append ((psym,count))
            psym = sym
            count = 1
        else:
            count += 1
    sc.append((psym,count))

    for sym, _ in sc:
        of.write(' {:3s}'.format(sym))
    of.write('\n')

    for _, count in sc:
        of.write(' {:3d}'.format(count))
    of.write('\n')

#writing the type of cell used
    if cartesian=='True':
    	of.write('Cartesian\n')
    else:
        of.write('Direct\n')
#writing the new coordinates
    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'
    count=0
    for iatom, atom in enumerate(new_coord):
        new_atom=dot(P_matrix,(atom-new_origin))
        print(new_atom)
        for dcoord in new_atom:
            of.write(cform % dcoord)
        of.write('\n')
        count += 1
    tot_natoms=count



atoms1,lat_scale,cell1X,cell1Y,cell1Z=read_POSCAR(str_name)

create_POSCAR(atoms1,lat_scale,cell1X,cell1Y,cell1Z,str_name,'False')



