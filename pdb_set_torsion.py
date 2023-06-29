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
if len(sys.argv)!=7:
	print('Usage pdb_set_torsion.py file.pdb atom1spec atom2spec atom3spec atom4spec angle')
	print('e.g. pdb_set_torsion.py bdglc.pdb 1.H5 1.C5 1.C6 1.O6 60.00')
	sys.exit(1)
	
PDBfilename=sys.argv[1]
atom1spec=sys.argv[2]
atom2spec=sys.argv[3]
atom3spec=sys.argv[4]
atom4spec=sys.argv[5]
settorsionangle=sys.argv[6]
print('REMARK Filename = ',PDBfilename)

# define blank dictionary as molecular data structure
mol1={}
# read PDB file and set elements within dictionary mol1
# outpput from function is number of atoms
pdb_natoms=read_pdb_file(mol1,PDBfilename)

## listmol(mol1)
# list molecule data
# print mol1
calc_atom_radius(mol1)
#listmol(mol1)
calcbonds(mol1)

atom1=atomspec_to_atom(mol1,atom1spec)
#print 'Atom selected =',atom1
if (atom1<=0):
	sys.exit(1)
atom2=atomspec_to_atom(mol1,atom2spec)
if (atom2<=0):
	sys.exit(1)
#print 'Atom selected =',atom2
atom3=atomspec_to_atom(mol1,atom3spec)
if (atom3<=0):
	sys.exit(1)
#print 'Atom selected =',atom3
atom4=atomspec_to_atom(mol1,atom4spec)
if (atom4<=0):
	sys.exit(1)
#print 'Atom selected =',atom4
print("REMARK Atom numbers for selected atoms =",atom1,atom2,atom3,atom4)

tors=torsion4atoms(mol1,atom1,atom2,atom3,atom4)
print('REMARK Initial torsion angle (degrees) =',tors)
#print "-------------------------------------------------"

# Now code for checking ring bond
# if not a ring bond then set torsion angle
#atom1name="1.C1"
#atom2name="1.O1"
#atom1number=atomspec_to_atom(mol1,atom1name)
#atom2number=atomspec_to_atom(mol1,atom2name)

#print "Testing for ring bond"
#listatom(mol1,atom2)
#listatom(mol1,atom3)
##print isbonded(mol1,2,3)

# NB Could check for 3 bonds
# atom1-atom2 atom2-atom3 atom3-atom4
# ring bond checker DOES check atom2-atom3 bonded

# walk graph to find atoms only on atom2 side of bond
# exit if ring if bond is part of a ring - can NOT set torsion
# NB ring_size_for_bond MUST be called before setting a torsion
# as it sets the atom attribute 'flag1' for atom z axis rotation or not
bondringsize=ring_size_for_bond(atom2,atom3,mol1)
if (bondringsize>0):
	print("Error: selected torsion is part of a ring system")
	print("Torsion angle cannot be set ")
	print("Ring zize for this bond =",bondringsize)
	listatom(mol1,atom2)
	listatom(mol1,atom3)
	sys.exit(1)
#print "Ring zize for this bond =",bondringsize
# if ringsize is zero then OK to perform rotation around z axis
# in order to achive desired torsion angle
# atoms OK to rotate have atom attribute 'flag1' greater than zero

# check which atoms have flag1 set and may be rotated as torsion is changed
# list_flag1(mol1)

# increase torsion to get to specified value (0 to 360 degrees)
newangle=float(settorsionangle)
oldangle=tors
change=newangle-oldangle
##print "CHANGE =",change
change_torsion(mol1,atom1,atom2,atom3,atom4,change)

# list all atom coordinates
# listmol(mol1)

# Now find new torsion angle
tors=torsion4atoms(mol1,atom1,atom2,atom3,atom4)
print('REMARK New torsion angle (degrees)=',tors)
#print "-------------------------------------------------"

##xvalue=3.0
##svalue="hello"
##print '%s %6.3f ' % (svalue,xvalue)

# write new coordinates to STDOUT
write_pdb_file(mol1)






