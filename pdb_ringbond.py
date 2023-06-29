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

# Program finds size of ring holding this bond
# NB If more than one such ring exists it is NOT guaranteed to find the smallest 

# command line arguments
# programname=sys.argv[0]
#print 'Number of arguments = ',len(sys.argv)
if len(sys.argv)!=4:
	print('Usage pdb_ringbond.py file.pdb atom1spec atom2spec')
	print('e.g. pdb_ringbond.py bdglc.pdb 1.C1 1.C2')
	sys.exit()
	
PDBfilename=sys.argv[1]
atom1spec=sys.argv[2]
atom2spec=sys.argv[3]
print('Filename = ',PDBfilename)
#print 'Atom 1 specification=',atom1spec
#print 'Atom 2 specification=',atom2spec
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
	sys.exit
atom2=atomspec_to_atom(mol1,atom2spec)
if (atom2<=0):
	sys.exit
#print 'Atom selected =',atom2
print("-------------------------------------------------")

# Now code for checking ring bond
# if not a ring bond then set torsion angle

atom1name="1.C1"
##atom2name="1.O1"
##atom1number=atomspec_to_atom(mol1,atom1name)
##atom2number=atomspec_to_atom(mol1,atom2name)
print("Testing for ring bond")
listatom(mol1,atom1)
listatom(mol1,atom2)
##print isbonded(mol1,2,3)
# walk graph to find atoms only on atom2 side of bond
bondringsize=ring_size_for_bond(atom1,atom2,mol1)
print("Ring size for this bond =",bondringsize)








