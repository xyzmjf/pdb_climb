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

# create a fragment file using specified atoms
# in the fragment file the join atom is (confusingly) first 
# and the deleted atom is second
# program needs to walk molecular graph (if needed) at delete atom end of bond
# flag up which atoms NOT to include in fragment
# then output join atom first, del-atom second and then the rest of the atoms

# command line arguments
# programname=sys.argv[0]
#print 'Number of arguments = ',len(sys.argv)
if len(sys.argv)!=4:
	print('Usage pdb_make_fragment.py file.pdb outer-atom inner-atom')
	print('e.g. pdb_make_fragment.py bdglc.pdb 1.O1 1.C1')
	sys.exit()
# NB inner atom is retained => first atom in output
# outer atom is retained => second atom in output
# but all atoms (if any) beyond the outer atom in the molecular graph are deleted

	
PDBfilename=sys.argv[1] 
# varous names for the same atom(s)
# delatom = outer-atom=atom1
# joinatom = inner-atom=atom2
delatom=sys.argv[2]
joinatom=sys.argv[3]

print('REMARK Filename = ',PDBfilename)
print('REMARK delatom=',delatom)
print('REMARK joinatom=',joinatom)


# define blank dictionary as molecular data structure
mol1={}
# read PDB file and set elements within dictionary mol1
# outpput from function is number of atoms
pdb_natoms1=read_pdb_file(mol1,PDBfilename)

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

#--------------------------------

#print "atom numbers=",atom1,atom2

# check atom1 and atom2 are bonded
if (not isbonded(mol1,atom1,atom2)):
	print("Error atoms not bonded",delatom,joinatom)
	print("Program must exit")
	sys.exit(1)

# check for ring containing bond atom1-atom2
# if bond in ring then program must exit
ringsize=ring_size_for_bond(atom1,atom2,mol1)
if (ringsize>0):
	print("Error ring bond found",joinatom,delatom)
	print("Program must exit")
	sys.exit(1)


#call to ring_size_for_bond sets 'flag1' attribute for atoms in one portion of molecular graph
reset_flag1(mol1)
ringsize=ring_size_for_bond(atom2,atom1,mol1)
##list_flag1(mol1)


# NB Not yet consistent atom numbers 1,2,3..... in output
write_pdb_atom(mol1,atom2)
write_pdb_atom(mol1,atom1)
##sys.exit(1)
# after walking molecular graph keep all further atoms with flag1 attribute set to zero
write_pdb_flag1_under(1,mol1)





