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
if len(sys.argv)!=3:
	print('Usage pdb_set_resname.py file.pdb name')
	print('e.g. pdb_set_resname.py bdglc.pdb GLC')
	sys.exit()
	
PDBfilename=sys.argv[1]
new_resname=sys.argv[2]
print('REMARK Filename = ',PDBfilename)

# define blank dictionary as molecular data structure
mol1={}
# read PDB file and set elements within dictionary mol1
# outpput from function is number of atoms
pdb_natoms=read_pdb_file(mol1,PDBfilename)

set_atom_resname(mol1,new_resname)
# write new coordinates to STDOUT
write_pdb_file(mol1)






