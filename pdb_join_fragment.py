#!/usr/bin/env python3


# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.

# THIS PROGRAM IS MADE AVAILABLE FOR DISTRIBUTION WITHOUT ANY FORM OF WARRANTY TO THE
# EXTENT PERMITTED BY APPLICABLE LAW. THE COPYRIGHT HOLDER PROVIDES THE PROGRAM AS IS
# WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM LIES
# WITH THE USER. SHOULD THE PROGRAM PROVE DEFECTIVE IN ANY WAY, THE USER ASSUMES THE
# COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION. THE COPYRIGHT HOLDER IS NOT
# RESPONSIBLE FOR ANY AMENDMENT, MODIFICATION OR OTHER ENHANCEMENT MADE TO THE PROGRAM
# BY ANY USER WHO REDISTRIBUTES THE PROGRAM SO AMENDED, MODIFIED OR ENHANCED.

# IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL THE
# COPYRIGHT HOLDER BE LIABLE TO ANY USER FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL,
# INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
# PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE
# OR LOSSES SUSTAINED BY THE USER OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO
# OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGES.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.




import sys
# access functions in mol_utils.py file
from mol_utils import *
from math import *




# ------------------start of main program----------------------

# command line arguments
# programname=sys.argv[0]
#print 'Number of arguments = ',len(sys.argv)
if len(sys.argv)!=5:
	print('Usage pdb_join_fragment.py file.pdb del-atom join-atom fragment.pdb')
	print('e.g. pdb_join_fragment.py bdglc.pdb 1.O1 1.C1 frag_4bdglc.pdb')
	sys.exit()
	
PDBfilename=sys.argv[1]
delatom=sys.argv[2]
joinatom=sys.argv[3]
FRAGPDBfilename=sys.argv[4]

print('REMARK Filename = ',PDBfilename)
print('REMARK delatom=',delatom)
print('REMARK joinatom=',joinatom)
print('REMARK Fragemnt PDB file=',FRAGPDBfilename)

# define blank dictionary as molecular data structure
mol1={}
mol2={}
# read PDB file and set elements within dictionary mol1
# outpput from function is number of atoms
pdb_natoms1=read_pdb_file(mol1,PDBfilename)
pdb_natoms2=read_pdb_file(mol2,FRAGPDBfilename)

# check that atoms exist
atom1=atomspec_to_atom(mol1,delatom)
#print 'Atom selected =',atom1
if (atom1<=0):
	print("Error: Atom not found",delatom)
	print("Program must exit")	
	sys.exit(1)
atom2=atomspec_to_atom(mol1,joinatom)
if (atom2<=0):
	print("Error: Atom not found",joinatom)
	print("Program must exit")	
	sys.exit(1)



#set_atom_resnumber(mol1,new_resnumber)
# write new coordinates to STDOUT
#write_pdb_file(mol1)

# list atom data
##listmol(mol1)
##listmol(mol2)

# compute bonds
calc_atom_radius(mol1)
calcbonds(mol1)
calc_atom_radius(mol2)
calcbonds(mol2)

#--------------------------------

#print "atom numbers=",atom1,atom2

# check atom1 and atom2 are bonded
if (not isbonded(mol1,atom1,atom2)):
	print("Error atoms not bonded",delatom,joinatom)
	print("Program must exit")
	sys.exit(1)

# check for ring containing bond atom1-atom2
# if bond in ring then program myst exit
ringsize=ring_size_for_bond(atom2,atom1,mol1)
if (ringsize>0):
	print("Error ring bond found",joinatom,delatom)
	print("Program must exit")
	sys.exit(1)

#------------------------
# do geometry stuff on mol1

# if not a ring bond then atom attribute 'flag1' is now > 1
# for atoms that must be removed from mol1
# list_flag1(mol1)
#list_flag1_under(2,mol1)

# translate mol1 so that join atom is at origin
new_origin_atom(mol1,atom2)

# find transform to place delatom (atom1) on z axis
trmat={}
##x1=mol1[atom1,'x']
##y1=mol1[atom1,'y']
##z1=mol1[atom1,'z']
findrot1atom(atom1,trmat,mol1)
##findrot1(x1,y1,z1,trmat)
# apply transformation to all atoms in mol1
rot1mol(mol1,trmat)

# list atom coordinates
#listmol(mol1)

#--------------------------------

# do geometry stuff on mol2
# translate mol2 so that join atom (first atom) is at origin
join2=1
new_origin_atom(mol2,join2)

# find transform to place mol2 deleted atom (second atom) on Z axis
trmat={}
del2=2
findrot1atom(del2,trmat,mol2)
##x1=mol2[del2,'x']
##y1=mol2[del2,'y']
##z1=mol2[del2,'z']
##findrot1(x1,y1,z1,trmat)

# apply transformation to all atoms in mol2
rot1mol(mol2,trmat)
xrotation(mol2,180.0)

# bond length now taken to be sum of radii for two atoms being joined
nameatom1=mol1[atom2,'atomname']
# molecule 2 assumed to be fragment so first atom is joining atom
nameatom2=mol2[1,'atomname']
bondradius1=bond_radius_from_atomname(nameatom1)
bondradius2=bond_radius_from_atomname(nameatom2)
bondlength=bondradius1+bondradius2

print("REMARK bond length for ",nameatom1," and ",nameatom2," set to ",bondlength)
#bondlength=1.5
ztranslation(mol2,bondlength)

# set 'flag1' attribute so that deleted atom (second atom) in fragment 
# is not included in combined molecule 
# value of flag1 < 2 implies atom is included in combined molecule
# value of flag1 < 1 implies atoms to be included in bumpcheck 
reset_flag1(mol2)
mol2[2,'flag1']=2
mol2[1,'flag1']=1


#=====

# list atom coordinates
##listmol(mol2)


##list_flag1(mol1)
##list_flag1(mol2)
# what are joining atoms in mol1 and mol2
# these need to be excluded from the bumpcheck
##a1=atom2
##a2=2
bumpdist=bump_check_two_mols_flag1_under(1,mol1,mol2)
print("REMARK Initial bump distance=",bumpdist)
best_angle=0

# Now rotate mol2 around z axis to maximise the minimum bump distance
for step in range (1,73):
	zangle=5.0
	zrotation(mol2,zangle)
	newbumpdist=bump_check_two_mols_flag1_under(1,mol1,mol2)
	##print "step and Shortest bump distance=",step,newbumpdist
	if (newbumpdist>bumpdist):
		bumpdist=newbumpdist
		best_angle=step*zangle
		##print "Best bump distance so far =",bumpdist
		##print "Best angle=",best_angle

# retain best coordinates - maximising minimum bump distance
zrotation(mol2,best_angle)
bumpdist=bump_check_two_mols_flag1_under(1,mol1,mol2)
print("REMARK Final bump distance=",bumpdist)
##sys.exit(1)

## CAN REMOVE LAST CHUNK OF CODE INTO A FUNCTION CALLED ZROTATION BUMP SCAN 

#--------------------------------
# now list PDB atom records for only those atoms in combined molecule


# NB Not yet consistent atom numbers 1,2,3..... in output
write_pdb_flag1_under(2,mol1)
write_pdb_flag1_under(2,mol2)



# PERL CODE FROM molbuild.pl below here
#=========================================================

# perform sanity checks 
# i.e. that two atoms form a bond, bond is not part of a ring
##sanity_check_bond($atom1mol1,$atom2mol1,\%mol1);
##sanity_check_bond($atom1mol2,$atom2mol2,\%mol2);

# print number of bonds for atoms forming new bond
##$nbonds1=nbonds($atom1mol1,\%mol1);
##$nbonds2=nbonds($atom1mol2,\%mol2);
##print "Numbers of bonds = $nbonds1 $nbonds2\n";

# some atoms in structure(s) must be deleted 
# for selected bond atom1-atom2 
# the atoms to be deleted are all atoms at or beyond atom2 
# side of this bond
# initially set 'deleted' flag to 0 for all atoms
##clear_deleted_flags(\%mol1);
##clear_deleted_flags(\%mol2);

# set deleted flags so that some atoms not present in output
##set_deleted_flag($atom2mol1,\%mol1);
##set_deleted_flag($atom2mol2,\%mol2);
# set deleted graph for all atoms beyond atom2 end of bond atom1-atom2
##walk_graph($atom1mol1,$atom2mol1,\%mol1);
##walk_graph($atom1mol2,$atom2mol2,\%mol2);



# list deleted flags
#list_deleted_flags(\%mol1);
#list_deleted_flags(\%mol2);

# translate molecule1 so that first selected atom at origin
##translate_to_origin($atom1mol1,\%mol1);
# compute rotation (transformation) matrix trmat
# print "DEBUG atom2mol1 = $atom2mol1 \n";
##compute_trmat($atom2mol1,\%trmat,\%mol1);
# rotate coordinates to place second selected atom on +z axis
##transform_coords(\%trmat,\%mol1);


# translate molecule2 so that first selected atom at origin
##translate_to_origin($atom1mol2,\%mol2);
# compute rotation (transformation) matrix trmat
##compute_trmat($atom2mol2,\%trmat,\%mol2);
# rotate coordinates to place second selected atom on +z axis
##transform_coords(\%trmat,\%mol2);
# rotate second selected atom onto -z axis
##xrotate(180.0,\%mol2);
# translate 2nd molecule along +z axis
##dist=(1.0*$bondlength);
##ztranslate($dist,\%mol2);



