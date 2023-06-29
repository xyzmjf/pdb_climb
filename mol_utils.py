

## utility functions for molecular structure code

import sys
from math import *


#------------------------------------------------------------
# function list_coords
# list atom cordinates for all atoms in molecule
#------------------------------------------------------------
def list_coords(mol):
	natoms=mol['natoms']
	print('Number of atoms = ',natoms)
	#print 'list of coordinates '
	for i in range(1,natoms+1):
		x=mol[i,'x']
		y=mol[i,'y']
		z=mol[i,'z']
		print('Atom index,x,y,z =',i,x,y,z)
# end of function

#------------------------------------------------------------
# function list_element_coords
# list element symbol and atom cordinates for all atoms in molecule
#------------------------------------------------------------
def list_element_coords(mol):
	natoms=mol['natoms']
	## print 'Number of atoms = ',natoms
	#print 'list of coordinates '
	for i in range(1,natoms+1):
		x=mol[i,'x']
		y=mol[i,'y']
		z=mol[i,'z']
		elementstr=mol[i,'element']
		##print elementstr,x,y,z
		print('%2s  %8.4f%8.4f%8.4f' % (elementstr,x,y,z))
# end of function



#------------------------------------------------------------
# function listmol
# list atom information for all atoms in molecule
#------------------------------------------------------------
def listmol(mol):
        natoms=mol['natoms']
        print('number of atoms = ',natoms)
        #print 'list of coordinates '
        for i in range(1,natoms+1):
                listatom(mol,i)
# end of function



#------------------------------------------------------------
# function listatom
# list coordinates, residue name, number etc for a single atom
#------------------------------------------------------------
def listatom(mol,i):
	xx=mol[i,'x']
	yy=mol[i,'y']
	zz=mol[i,'z']
	rname=mol[i,'resname']
	rnumber=mol[i,'resnumber']
	aname=mol[i,'atomname']
	print('Atom ',i,aname, rname,rnumber,xx,yy,zz)
# end of function




#------------------------------------------------------------
# function distsqatoms
# compute distance squared between 2 atoms in molecule
#------------------------------------------------------------
def distsqatoms(mol,i,j):
	x1=mol[i,'x']
	y1=mol[i,'y']
	z1=mol[i,'z']
	x2=mol[j,'x']
	y2=mol[j,'y']
	z2=mol[j,'z']
	distsq=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
	return(distsq)	
# end of function



#------------------------------------------------------------
# function distatoms
# compute distance between 2 atoms in molecule
# calls function distsqatoms
#------------------------------------------------------------
def distatoms(mol,i,j):
	distsq=distsqatoms(mol,i,j)
	dist=sqrt(distsq)
	return(dist)
# end of function


#------------------------------------------------------------
# function makebond
# define  bonding elements of mol data structure
# build up data for both i-j and j-i bonds
#------------------------------------------------------------
def makebond(mol,i,j):
	# current numbers of bonds for atoms i and j
	nbondsi=mol[i,'nbonds']
	nbondsj=mol[j,'nbonds']
	nbondedi=1+nbondsi
	nbondedj=1+nbondsj	
	mol[i,'nbonds']=nbondedi
	mol[j,'nbonds']=nbondedj
	mol[i,'bonded',nbondedi]=j
	mol[j,'bonded',nbondedj]=i
# end of function


#------------------------------------------------------------
# function isbonded
# is atom2 bonded to atom1 
# return 0 or 1
# by checking all atoms bonded to atom1
# NOT doing reverse check
#------------------------------------------------------------
def isbonded(mol,atom1,atom2):
	value=0
	nbonds1=mol[atom1,'nbonds']
	for j in range (1,nbonds1+1):
		atomj=mol[atom1,'bonded',j]
		if (atomj==atom2):
			value=1
	return(value)
# end of function


#------------------------------------------------------------
# function reset_flag1
# set flag1 to zero for all atoms in molecule
#------------------------------------------------------------
def reset_flag1(mol):
# list for all atoms
		natoms=mol['natoms']
		for i in range (1,natoms+1):
			mol[i,'flag1']=0
#			print "Atom and flag ",i,mol[i,'flag1']


#------------------------------------------------------------
# function list_flag1
# list flag1 values for all atoms in molecule
#------------------------------------------------------------
def list_flag1(mol):
# list for all atoms
		natoms=mol['natoms']
		for i in range (1,natoms+1):
			print("------")
			print("Atom and flag ",i,mol[i,'flag1'])
			listatom(mol,i)


#------------------------------------------------------------
# function list_flag1_under
# set atoms where attribute flag1 < value
#------------------------------------------------------------
def list_flag1_under(value,mol):
# list all atoms where flag1 < value
		natoms=mol['natoms']
		for i in range (1,natoms+1):
			if (mol[i,'flag1']<value):
				print("------")
				print("Atom and flag ",i,mol[i,'flag1'])
				listatom(mol,i)


#------------------------------------------------------------
# function write_pdb_flag1_under
# write PDB record for atoms where attribute flag1 < value
#------------------------------------------------------------
def write_pdb_flag1_under(value,mol):
# list all atoms where flag1 < value
		natoms=mol['natoms']
		for i in range (1,natoms+1):
			if (mol[i,'flag1']<value):
				##print "------"
				##print "Atom and flag ",i,mol[i,'flag1']
				##listatom(mol,i)
				write_pdb_atom(mol,i)




#------------------------------------------------------------
# function walk_mol_graph
# walk molecule graph at atom2 end of bond atom1-atom2
# set flag1 for atoms on this side of graph
# later use for rotation deletion etc.
# NB Recursion depth issue if you give a ring bond....
# NB Ensure reset_flag1 called before calling this function
#------------------------------------------------------------
##def walk_mol_graph(mol,atom1,atom2,depth):
##	if isbonded(mol,atom1,atom2):
# now walk graph on atom2 side
##		nbonds2=mol[atom2,'nbonds']
##		if (nbonds2>1):	
##			for j in range (1,nbonds2+1):
# find atom number of atoms bonded to atom2
# if flag not already set
# then set flag and call function 
# recursively for bonded atoms
##				atomj=mol[atom2,'bonded',j]
##				if (atomj<>atom1):
##					if (mol[atomj,'flag1']==0):
##						newdepth=1+depth
##						listatom(mol,atomj)
##						print "Depth =  ",newdepth
##						mol[atomj,'flag1']=1
##						print "Flag  ",mol[atomj,'flag1']
##						walk_mol_graph(mol,atom2,atomj,newdepth)
##	else:
##		print "ERROR: atoms not bonded ",atom1,atom2
# end function

#------------------------------------------------------------
# function is_ring_bond
# recursively walk molecule graph at atom2 end of bond atom1-atom2
# if depth>2 and we get back to atom1 then must be part of a ring
# Initially call function with depth set=1
# NB NOT CURRENTLY USED
#------------------------------------------------------------
##def is_ring_bond(mol,atom1,atom2,recdepth):
# recdepth is depth if recursion
# set flag1 to zero for all atoms in molecule
# then set flag1 = -1 for atom1
# and set flag1 = +1 for atom2
##	if (recdepth==0):
##		reset_flag1(mol)
##		mol[atom1,'flag1']=(-1)
##		mol[atom2,'flag1']=1
##	atom1flag=mol[atom1,'flag1']
##	if (atom1flag>0):
##		print "--------------------------------"
##		print "Found first atom again"
##		print "Depth=",recdepth
##	print "Recursion depth =",recdepth
# find number of bonds for atom2 
##	nbonds2=mol[atom2,'nbonds']
##	if (nbonds2>=1):
# if atom2 has 4 bonds then j =1,2,3,4	
##		for j in range (1,nbonds2+1):
# find atoms bonded to atom2
# if flag not already set
# then set flag and call function 
# recursively for bonded atoms
# atomj is atom index number for atom bonded to atom2
##			atomj=mol[atom2,'bonded',j]
##			flagj=mol[atomj,'flag1']
# if flag not set for this atom then set it
# increase depth and call function recursively
##			if (flagj==0):
##				newdepth=1+recdepth
##				listatom(mol,atomj)
##				print "Depth =  ",newdepth
##				mol[atomj,'flag1']=newdepth	
##				is_ring_bond(mol,atom1,atomj,newdepth)
##	print "Atom1flag=",atom1flag
# end function			


#------------------------------------------------------------
# function ring_size
# just sets initial parameters and calls function check_ring_bond
# ring size stored in molecule attribute 'ringsize'
# ------------------------------------------------------------
def ring_size_for_bond(atom1,atom2,mol):
# set molecule attribute 'ringsize' to zero
# later set in ring bond checking code
# used to indicate size of ring found for this bond
	mol['ringsize']=0
# set atom attribute 'flag1' to zero for all atoms
# then set to 1 and 2 for bond atoms
# this attribute is used to track bond separation depth in recursive code
# Do not check if initial atoms not bonded..
	if (isbonded(mol,atom1,atom2)):
		reset_flag1(mol)
		mol[atom1,'flag1']=1
		mol[atom2,'flag1']=2
		check_ring_bond(atom1,atom2,mol)
	else:
		print("Error: Atoms not bonded ",atom1,atom2)
		print("Program will exit !")
		exit(0)
	return(mol['ringsize'])

#------------------------------------------------------------
# function check_ring_bond
# walk along connectivity graph
# check all atoms on bondatom2 side of bondatom1-bondatom2 bond
# if cycles round and back to atom $atm1 then starting bond
# must be part of a ring 
# to start check call with atoms atom1 atom1 atom2
# ring size stored in molecule attribute 'ringsize'
# ------------------------------------------------------------
def check_ring_bond(atom1,atom2,mol):
	nbonds2=mol[atom2,'nbonds']
	# $depth2 was set to 2 before recursion
	depth2=mol[atom2,'flag1']
# only check if atom 2 has more than 1 bond
# If only 1 bond then no further bonds at far end
	if (nbonds2>1):
		##print "Checking atom "
		##listatom(mol,bondatom2)
		##print "Depth=",depth2
		#listatom(mol,atm1)
		for bond in range (1,nbonds2+1):
		# if atom2 has 4 bonds then j=1,2,3,4
			j=mol[atom2,'bonded',bond]
			depthj=mol[j,'flag1']
			##namej=$$refmol{$j,'atomname'};
			#
			# if starting atom atom1 is found
			# and depth of search is greater than 2 bonds
			# then bond chosen must be part of a ring
			#
			if (j==atom1 and depth2>2):
				##print "ALERT found bond "
				##print "in ring structure "
				##print "Ring size = ",depth2
				mol['ringsize']=depth2
			# depthj==0 ensures no backtracking allowed 
			if (depthj==0):
				# increase path depth for newly found atoms
				newdepth=1+depth2
				##print "Atom $j $namej depth=$newdepth \n";
				mol[j,'flag1']=newdepth;
				##print "Now check bond $bondatom2 $j\n";
				check_ring_bond(atom1,j,mol)
# end def





#------------------------------------------------------------
# function centremol
# Translate all atoms by subtracting geometric centre
# ------------------------------------------------------------
def centremol(mol):
	natoms=mol['natoms']
	xsum=0.0
	ysum=0.0
	zsum=0.0
	for i in range (1,natoms+1):
		xx=mol[i,'x']
		yy=mol[i,'y']
		zz=mol[i,'z']
		xsum=xsum+xx
		ysum=ysum+yy
		zsum=zsum+zz
	xmean=xsum/(float(natoms))
	ymean=ysum/(float(natoms))
	zmean=zsum/(float(natoms))
	
	# now subtract mean coordinates
	for i in range (1,natoms+1):
		xx=mol[i,'x']
		yy=mol[i,'y']
		zz=mol[i,'z']
		mol[i,'x']=xx-xmean
		mol[i,'y']=yy-ymean
		mol[i,'z']=zz-zmean
# end of function


#------------------------------------------------------------
# function listbonds
# List all bonds in molecule
# calls function listbonds1atom
# ------------------------------------------------------------
def listbonds(mol):
	natoms=mol['natoms']
	for i in range (1,natoms+1):
		listbonds1atom(mol,i)
# end of function


#------------------------------------------------------------
# function listbonds1atom
# List bonds for 1 atom in molecule
# ------------------------------------------------------------
def listbonds1atom(mol,i):
	nbondsi=mol[i,'nbonds']
	print() 
	listatom(mol,i)
	print('Number of bonded atoms =',nbondsi)	
	for j in range (1,nbondsi+1):
		atomj=mol[i,'bonded',j]
		listatom(mol,atomj)
		print('distance=',distatoms(mol,i,atomj))
# end of function




#------------------------------------------------------------
# function calcbonds
# Use atom radii to compute and save bonds in molecule
# ------------------------------------------------------------
def calcbonds(mol):
# slow calculation of likely bonds
	#cutoff=1.75
	natoms=mol['natoms']
# set initial numbers of bonds to zero
	for i in range(1,natoms+1):
		mol[i,'nbonds']=0
# now calculate bonds based on distance criteria
	for i in range(1,1+natoms):
		for j in range (1,1+natoms):
			if (j>i):			
				distsq=distsqatoms(mol,i,j)
				cutoff=mol[i,'radius']+mol[j,'radius']
				if (distsq<(cutoff*cutoff)):
# now populate bonding elements of data structure
# make i-j and j-i bond
					makebond(mol,i,j)
# end of function



#------------------------------------------------------------
# function calc_atom_radius
# simple name based calculation of atom bonding radius
# ------------------------------------------------------------
def calc_atom_radius(mol):
	# set atom radii for all atoms in molecule
	natoms=mol['natoms']
	# set initial colours, will hold for non identified atoms
	for i in range(1,1+natoms):
		mol[i,'radius']=0.8
	# now set radii for specific atoms
	for i in range(1,1+natoms):
		name=mol[i,'atomname']
		#print len(name)
# need to deal with case of protons named
# 1H 2H etc. 
		name0=name[0]
# evaluating name[1] i.e. 2nd character in name
# throws an error is name is only of length 1
		#name1=name[1]		
		if name0=='H':
			mol[i,'radius']=0.5
		if name0=='C':
			mol[i,'radius']=0.85
		if name0=='O':
			mol[i,'radius']=0.8
		if name0=='N':
			mol[i,'radius']=0.8
		if name0=='S':
			mol[i,'radius']=1.0
		if name0=='P':
			mol[i,'radius']=1.0
		#print 'first two characters in atom name = ',name0,name1
		#print 'Radius=',mol[i,'radius']






#------------------------------------------------------------
# function bond_radius_from_atomname
# simple name based calculation of atom bonding radius
# ------------------------------------------------------------
def bond_radius_from_atomname(atomname):
	# set bonding radius used to calculate new bond lengths
	# Assuming single bonds at present
	# See table on Wikipedia http://en.wikipedia.org/wiki/Covalent_radius
	# 
		name=atomname
		
# need to deal with case of protons named
# 1H 2H etc. 

		# print name[0:2]
		name0=name[0]
		if (len(name)>=2):	
			name1=name[1]
#		default value rather large for some, small for others
		value=1.5

# evaluating name[1] i.e. 2nd character in name
# throws an error is name is only of length 1		
		if name0=='H':
			value=0.31
		# Check for names like 1H 1HA 1H2 etc.
		if (len(name)>=2):
			if (name[0].isdigit() and name[1]=='H'):
				value=0.31			
		if name0=='C':
			value=0.76
		if name0=='B':
			value=0.84
		if name0=='O':
			value=0.66
		if name0=='N':
			value=0.71
		if name0=='F':
			value=0.57
		if name0=='S':
			value=1.05
		if name0=='P':
			value=1.07
		if (name[0:2]=='Cl'):
			value=1.02
		if (name[0:2]=='Cl'):
			value=1.20
		if (name[0:2]=='Br'):
			value=1.20
		#print 'first two characters in atom name = ',name0,name1
		#print 'Radius=',mol[i,'radius']
		return(value)



#------------------------------------------------------------
# function element_from_atomname
# simple name based calculation element string
# ------------------------------------------------------------
def element_from_atomname(atomname):
		name=atomname
# need to deal with case of protons named
# 1H 2H etc. 
		# print name[0:2]
		## print "name and length=",name,len(name)
		name0=name[0]
		if (len(name)>=2):
			name1=name[1]
#		set an initial value, throw an error near exit if still set 
		value='XX'
# evaluating name[1] i.e. 2nd character in name
# throws an error is name is only of length 1		
		if name0=='H':
			value='H'
		# Check for names like 1H 1HA 1H2 etc.
		if (len(name)>=2):
			if (name[0].isdigit() and name[1]=='H'):
				value='H'			
		if name0=='C':
			value='C'
		if name0=='B':
			value='B'
		if name0=='O':
			value='O'
		if name0=='N':
			value='N'
		if name0=='F':
			value='F'
		if name0=='S':
			value='S'
		if name0=='P':
			value='P'
		if (name[0:2]=='Cl'):
			value='Cl'
		if (name[0:2]=='Br'):
			value='Br'
		#print 'first two characters in atom name = ',name0,name1
		#print 'Radius=',mol[i,'radius']
		if (value=='XX'):
			print("ERROR: Element symbol not set for atom ",atomname)
			print("Program must exit")
			sys.exit(1)
		return(value)





#------------------------------------------------------------
# function atomspec_to_atom
# find atom based on text specification
# example 113.CA for residue number 113 atom called CA
# ------------------------------------------------------------
def atomspec_to_atom(mol,spec):
	# split into residue number spec and atom name spec 
	# accepting period character as delimiter
	# check that the spec contains a period character
	found=spec.find('.')
	#print found
	if (found<0):
		print('ERROR: Character . needed in atom specifications',spec)
		sys.exit
	atomfound=0
	#print 'Specification =',spec
	natoms=mol['natoms']
	length=len(spec)
	resnumspec=spec[0:found]
	atomspec=spec[1+found:length]
	#[resnumspec,atomspec]=spec.split('.')
	#print 'Number of atoms=',natoms
	for i in range(1,1+natoms):		
		name1=mol[i,'atomname']
		resnum1=mol[i,'resnumber']
		#print 'Atoms ',name1,spec,len(name1),len(spec)
		if (name1==atomspec):			
			if (resnum1==resnumspec):
				atomfound=i
				#print 'Found atom with spec=',spec
				#print 'Atom number =',i
	if atomfound==0:
		print('Could not find atom with specification ',spec)
		sys.exit
	return(atomfound)
# end of function



#------------------------------------------------------------
# function read_pdb_file
# Reads ATOM and HETATM records from PDB format file
# stores relevant data in dictionary named mol
# most work done in function parse_pdb_line 
# ------------------------------------------------------------
def read_pdb_file(mol,filename):
#print 'If you get this far then should be OK to open file ',PDBfilename
	pdbfile=open(filename,'r')
	pdb_patoms=0
	pdb_hetatoms=0	
	atomcount=0
	for line in pdbfile.readlines():
	# print line
		line4=line[:4]
		line6=line[:6]
	# print line4,line6
		if line4=='ATOM':
			atomcount=atomcount+1
			pdb_patoms=pdb_patoms+1
	# call function to atom data from line
			parse_pdb_line(mol,atomcount,line)
		elif line6=='HETATM': 
			atomcount=atomcount+1
			pdb_hetatoms=pdb_hetatoms+1
	# call function to atom data from line
			parse_pdb_line(mol,atomcount,line)
	mol['natoms']=pdb_patoms+pdb_hetatoms
#print 'Number of ATOM records =',pdb_patoms
#print 'Number of HETATM records=',pdb_hetatoms
	pdb_natoms=pdb_patoms+pdb_hetatoms
	return(pdb_natoms)



#------------------------------------------------------------
# function findrot1 
# calculate transformation matrix that puts coordinate xx1,yy1,zz1
# onto positive z axis  
# ------------------------------------------------------------
def findrot1(xx1,yy1,zz1,trmat):
# Need to check 
# find sine and cosine of rotation angles
# needed to place point in zy plane (where x equals zero)
# if point already on z axis xx=yy=0 then old code had 
# divide by zero error
# check if xx squared is zero, if so then already in y-z plane
# so that rotation angle of zero required
# then c1=cos(0)=1 and s1=sin(0)=0

	if (xx1*xx1>0):
		c1=yy1/sqrt(xx1*xx1+yy1*yy1)
		s1=xx1/sqrt(xx1*xx1+yy1*yy1)
	else:
		c1=1
		s1=0
	ynew=s1*xx1+c1*yy1
	s2=ynew/sqrt(ynew*ynew+zz1*zz1)
	c2=zz1/sqrt(ynew*ynew+zz1*zz1)
	trmat[1,1]=c1
	trmat[1,2]=(-1.0)*s1
	trmat[1,3]=0.0
	trmat[2,1]=c2*s1
	trmat[2,2]=c1*c2
	trmat[2,3]=(-1.0)*s2
	trmat[3,1]=s1*s2
	trmat[3,2]=c1*s2
	trmat[3,3]=c2


#------------------------------------------------------------
# function findrot1atom 
# calculate transformation matrix that puts coordinate selected atom
# onto positive z axis  
# ------------------------------------------------------------
def findrot1atom(atom1,trmat,mol):
# Need to check 
# find sine and cosine of rotation angles
# needed to place point in zy plane (where x equals zero)
# if point already on z axis xx=yy=0 then old code had 
# divide by zero error
# check if xx squared is zero, if so then already in y-z plane
# so that rotation angle of zero required
# then c1=cos(0)=1 and s1=sin(0)=0
	xx1=mol[atom1,'x']
	yy1=mol[atom1,'y']
	zz1=mol[atom1,'z']
	findrot1(xx1,yy1,zz1,trmat)




#------------------------------------------------------------
# function torsion2_0360
# now selected atoms 2 and 3 are aligned along z axis 
# and x,y coords for sected atoms 1 and 4 are are xx1,yy1,xx2,yy2
# find torsion angle between two vectors projected
# onto the z axis
# Return torsion in range 0..360 degrees
#------------------------------------------------------------

def torsion2_0360(xx1,yy1,xx2,yy2):
	ra=sqrt(xx1*xx1+yy1*yy1)
	c1=yy1/ra
	s1=xx1/ra
# NB CONVENTION used in function torsion2_0360
# MUST MATCH THAT IN Z ROTATION FUNCTIONS 
# SUCH AS zrot1atom
# i.e. x'=x cos theta + y sin theta ; y'=y cos theta - x sin theta 
# OR the opposite convention
# x ' = x cos theta - y sin theta ; y'=ycos theta + x sin theta
# Currently applying the latter
	xnew=xx2*c1-yy2*s1
	ynew=xx2*s1+yy2*c1
	rd=sqrt(xnew*xnew+ynew*ynew)
	theta=(180.0/3.1415926)*acos(ynew/rd)
# now ensure in range 0 to 360 
	if (xnew<0.0): 
		theta=(-1.0)*theta
	theta=(-1.0)*theta
	if (theta<0.0):
		diff180=180.0-fabs(theta)
		theta=180.0+diff180
	return(theta)



#------------------------------------------------------------
# function rot1mol
# carry out rotation defined by matrix trmat
# on all atoms in molecule
#------------------------------------------------------------

def rot1mol(mol,trmat):
	natoms=mol['natoms']
	for i in range(1,natoms+1):
		xold=mol[i,'x']
		yold=mol[i,'y']
		zold=mol[i,'z']	
		xnew=trmat[1,1]*xold+trmat[1,2]*yold+trmat[1,3]*zold
		ynew=trmat[2,1]*xold+trmat[2,2]*yold+trmat[2,3]*zold
		znew=trmat[3,1]*xold+trmat[3,2]*yold+trmat[3,3]*zold
		mol[i,'x']=xnew
		mol[i,'y']=ynew
		mol[i,'z']=znew
	



#------------------------------------------------------------
# function rot1
# carry out rotation defined by matrix trmat
# on arrays xx,yy,zz
#------------------------------------------------------------

def rot1(i,xx,yy,zz,trmat):
	xnew=trmat[1,1]*xx[i]+trmat[1,2]*yy[i]+trmat[1,3]*zz[i]
	ynew=trmat[2,1]*xx[i]+trmat[2,2]*yy[i]+trmat[2,3]*zz[i]
	znew=trmat[3,1]*xx[i]+trmat[3,2]*yy[i]+trmat[3,3]*zz[i]
	xx[i]=xnew
	yy[i]=ynew
	zz[i]=znew





#------------------------------------------------------------
# function torsion4atoms
# return torsion angle for four atoms at1,at2,at3,at4
# ------------------------------------------------------------
def torsion4atoms(mol,at1,at2,at3,at4):
	i1=1
	i2=2
	i3=3
	x1=mol[at1,'x']
	y1=mol[at1,'y']
	z1=mol[at1,'z']
	x2=mol[at2,'x']
	y2=mol[at2,'y']
	z2=mol[at2,'z']	
	x3=mol[at3,'x']
	y3=mol[at3,'y']
	z3=mol[at3,'z']
	x4=mol[at4,'x']
	y4=mol[at4,'y']
	z4=mol[at4,'z']
	xx={}
	yy={}
	zz={}
	trmat={}
	xx[1]=x1-x2
	yy[1]=y1-y2
	zz[1]=z1-z2
	xx[2]=x3-x2
	yy[2]=y3-y2
	zz[2]=z3-z2
	xx[3]=x4-x2
	yy[3]=y4-y2
	zz[3]=z4-z2
	xx2=xx[2]
	yy2=yy[2]
	zz2=zz[2]
# find rotation that places atom 3 on z axis
	findrot1(xx2,yy2,zz2,trmat)
	#print trmat
#   carry out this rotation on atoms 1 and 4
#   hence elements 1,3 of arrays xx,yy,zz
	rot1(i1,xx,yy,zz,trmat)
	rot1(i2,xx,yy,zz,trmat)
	#print 'i3=',i3
	#print 'xx=',xx
	#print 'yy=',yy
	#print 'zz=',zz
	rot1(i3,xx,yy,zz,trmat)
	angle=torsion2_0360(xx[1],yy[1],xx[3],yy[3])
	#print 'torsion angle=',angle	
	return(angle)



#------------------------------------------------------------
# function change_torsion
# Alter torsion angle for four atoms at1,at2,at3,at4
# NB Unlike the torsion4atoms measuring function
# this once actually modifies coordinates - of all atoms
# actions => translate molecule such that at2 at origin
# 	compute rotation to place atom at3 on z axis
# 	apply rotation to all atoms
# 	compute torsion as angle vetween at1-at2 and at3-at4 bonds
#	when viewed along z axis
#	Then perform extra z axis rotation, but only for atoms
#	where atom attribute 'flag1' greater than zero
# The 'flag1' attribute is set by calling the function
# ring_size_for_bond
# ------------------------------------------------------------
def change_torsion(mol,at1,at2,at3,at4,changeangle):
	natoms=mol['natoms']
	new_origin_atom(mol,at2)
# NB EMPIRICAL !!!!
	# changeangle=(-1.0)*changeangle
# WAS needed when x,y conventions upon z-rotation were NOT consistent
	# print "In change_torsion function changeangle=",changeangle
# translate all atoms so that atom at2 at origin
	##x2=mol[at2,'x']
	##y2=mol[at2,'y']
	##z2=mol[at2,'z']	
	##for i in range(1,natoms+1):
##		xnew=mol[i,'x']-x2
##		ynew=mol[i,'y']-y2
##		znew=mol[i,'z']-z2
##		mol[i,'x']=xnew
##		mol[i,'y']=ynew
##		mol[i,'z']=znew
# find transformation to put atom3 on + z axis
	x3=mol[at3,'x']
	y3=mol[at3,'y']
	z3=mol[at3,'z']
	trmat={}
	xx={}
	yy={}
	zz={}
	findrot1(x3,y3,z3,trmat)
# apply transformation to all atoms
# Note that function rot1 only operates on array elements
# TRY USING rot1mol !!!!!!!
	for i in range(1,natoms+1):
		xx[1]=mol[i,'x']
		yy[1]=mol[i,'y']
		zz[1]=mol[i,'z']
		rot1(1,xx,yy,zz,trmat)
		mol[i,'x']=xx[1]
		mol[i,'y']=yy[1]
		mol[i,'z']=zz[1]
# now find transformed x,y coordinates of atom1 and atom4
	xx[1]=mol[at1,'x']
	yy[1]=mol[at1,'y']
	xx[2]=mol[at4,'x']
	yy[2]=mol[at4,'y']
	angle=torsion2_0360(xx[1],yy[1],xx[2],yy[2])
	##print 'change in torsion angle=',changeangle	
	#rotationangle=30.0
	##print "Coordinates before zrotation_flagged_atoms"
	##listmol(mol)	
	zrotation_flagged_atoms(mol,changeangle)


#------------------------------------------------------------
# function zrotation_flagged_atoms	
# perform z axis rotation for atoms where
# 'flag1' attribute > 0
# ------------------------------------------------------------
def zrotation_flagged_atoms(mol,angle):
	##print "Angle of rotation=",angle
	radangle=angle*3.1415926/180.0
	## print "Rad angle=",radangle	
	s=sin(radangle)
	c=cos(radangle)
	## print "sin cos =",s,c
	natoms=mol['natoms']
	for i in range(1,natoms+1):
		flagi=mol[i,'flag1']
		if (flagi>0):
			#print "Z axis rotation for atom ",i
			##listatom(mol,i)
			zrot1atom(i,s,c,mol)
			##listatom(mol,i)


#------------------------------------------------------------
# functions xrot1atom
# rotate a single atom about x axis
# s,c are sin and cos of rotation angle
# ------------------------------------------------------------
def xrot1atom(i,s,c,mol):
# positive x axis rotation of atom i
#   s,c are sin and cos of rotation angle
	zi=mol[i,'z']
	yi=mol[i,'y']
	znew=zi*c-yi*s
	ynew=yi*c+zi*s
	mol[i,'y']=ynew
	mol[i,'z']=znew


#------------------------------------------------------------
# functions yrot1atom
# rotate a single atom about y axis
# s,c are sin and cos of rotation angle
# ------------------------------------------------------------
def yrot1atom(i,s,c,mol):
# positive x axis rotation of atom i
#   s,c are sin and cos of rotation angle
	xi=mol[i,'x']
	zi=mol[i,'z']
	xnew=xi*c+zi*s
	znew=zi*c-xi*s
	mol[i,'x']=xnew
	mol[i,'z']=znew



#------------------------------------------------------------
# functions zrot1atom
# function used by bond_rotation function
# they rotate a single atom about z axis
# s,c are sin and cos of rotation angle
# ------------------------------------------------------------
def zrot1atom(i,s,c,mol):
# positive z axis rotation of atom i
#   s,c are sin and cos of rotation angle
# NB CONVENTION used in function torsion2_0360
# MUST MATCH THAT IN Z ROTATION FUNCTIONS 
# SUCH AS zrot1atom
# i.e. x'=x cos theta + y sin theta ; y'=y cos theta - x sin theta 
# OR the opposite convention
# x ' = x cos theta - y sin theta ; y'=ycos theta + x sin theta
# Currently applying the latter
	xi=mol[i,'x']
	yi=mol[i,'y']
	xnew=xi*c-yi*s
	ynew=yi*c+xi*s
	mol[i,'x']=xnew
	mol[i,'y']=ynew
	## print "in z rot1 atom",i,xnew,ynew




#------------------------------------------------------------
# functions xrotation
# rotate all atoms in molecule
# they rotate about x axis
# s,c are sin and cos of rotation angle
# ------------------------------------------------------------
def xrotation(mol,angle):
	#print "Angle of rotation=",angle
	radangle=angle*3.1415926/180.0
	s=sin(radangle)
	c=cos(radangle)
	#print "sin cos =",s,c
	natoms=mol['natoms']
	for i in range(1,natoms+1):
		#print "x axis rotation for atom ",i
		#listatom(mol,i)
		xrot1atom(i,s,c,mol)
		#listatom(mol,i)



#------------------------------------------------------------
# functions yrotation
# rotate all atoms in molecule
# they rotate about y axis
# s,c are sin and cos of rotation angle
# ------------------------------------------------------------
def yrotation(mol,angle):
	#print "Angle of rotation=",angle
	radangle=angle*3.1415926/180.0
	s=sin(radangle)
	c=cos(radangle)
	#print "sin cos =",s,c
	natoms=mol['natoms']
	for i in range(1,natoms+1):
		#print "y axis rotation for atom ",i
		#listatom(mol,i)
		yrot1atom(i,s,c,mol)
		#listatom(mol,i)


#------------------------------------------------------------
# functions zrotation
# rotate all atoms in molecule
# they rotate about z axis
# s,c are sin and cos of rotation angle
# ------------------------------------------------------------
def zrotation(mol,angle):
	#print "Angle of rotation=",angle
	radangle=angle*3.1415926/180.0
	s=sin(radangle)
	c=cos(radangle)
	#print "sin cos =",s,c
	natoms=mol['natoms']
	for i in range(1,natoms+1):
		#print "Z axis rotation for atom ",i
		#listatom(mol,i)
		zrot1atom(i,s,c,mol)
		#listatom(mol,i)






#------------------------------------------------------------
# function parse_pdb_line
# find coordinates, residue name and number for single ATOM or HETATM record in PDB file
# ------------------------------------------------------------
# define a function - do not forget the colon:
def parse_pdb_line(mol,atomnumber,line):
	#print 'In function parse_pdb_line'
	xstr=line[30:38]
	ystr=line[39:47]
	zstr=line[48:56]
	atmname=line[12:17]
	atmname2=atmname.strip()	
	#print atmname,atmname2,xstr, ystr, zstr
	mol[atomnumber,'x']=float(xstr)
	mol[atomnumber,'y']=float(ystr)
	mol[atomnumber,'z']=float(zstr)
	mol[atomnumber,'atomname']=atmname2
	mol[atomnumber,'resname']=line[17:20].strip()
	mol[atomnumber,'resnumber']=line[24:27].strip()
# end of function


#------------------------------------------------------------
# function write_pdb_file
# writes to STDOUT
# write all ATOM/HETATM records in PDB file format
# ------------------------------------------------------------
def write_pdb_file(mol):
	natoms=mol['natoms']
	for i in range(1,1+natoms):
		write_pdb_atom(mol,i)


#------------------------------------------------------------
# function write_pdb_atom
# write single ATOM/HETATM record in PDB file format
# NB Format string lifted from mdxvu C code
# NB NOT currently reading chain ID
# NOT dealing with fields after z coordinate
# ------------------------------------------------------------
def write_pdb_atom(mol,atom):
	xcoord=mol[atom,'x']
	ycoord=mol[atom,'y']
	zcoord=mol[atom,'z']
	resname=mol[atom,'resname']
	resnumber=mol[atom,'resnumber']
	atomname=mol[atom,'atomname']
	chainid=" "
	print('ATOM  %5d  %-4.4s%3.3s %c %3.3s    %8.3f%8.3f%8.3f' % (atom,atomname,resname,chainid,resnumber,xcoord,ycoord,zcoord))


#------------------------------------------------------------
# function set_atom_resname
# set residue name to single value for ALL atoms in molecule
# may later write PDB format file for this single residue
# ------------------------------------------------------------
def set_atom_resname(mol,name):
	natoms=mol['natoms']
	for i in range(1,1+natoms):
		mol[i,'resname']=name


#------------------------------------------------------------
# function set_atom_resnumber
# set residue numbers to single value for ALL atoms in molecule
# may later write PDB format file for this single residue
# ------------------------------------------------------------
def set_atom_resnumber(mol,name):
	natoms=mol['natoms']
	for i in range(1,1+natoms):
		mol[i,'resnumber']=name



#------------------------------------------------------------
# function new_origin_coordinate
# translate ALL atoms in molecule
# so that new coordinate origin is at x0,y0,z0
# ----------------------------------------------------------
def new_origin_coordinate(mol,x0,y0,z0):
	natoms=mol['natoms']
# translate all atoms so that x0,y0,z0 at origin	
	for i in range(1,natoms+1):
		xnew=mol[i,'x']-x0
		ynew=mol[i,'y']-y0
		znew=mol[i,'z']-z0
		mol[i,'x']=xnew
		mol[i,'y']=ynew
		mol[i,'z']=znew
	
#------------------------------------------------------------
# function new_origin_atom
# translate ALL atoms in molecule
# so that new coordinate origin is centred on atom1
# ----------------------------------------------------------
def new_origin_atom(mol,atom1):
	natoms=mol['natoms']
	x0=mol[atom1,'x']
	y0=mol[atom1,'y']
	z0=mol[atom1,'z']
# translate all atoms so that atom1 at origin	
	for i in range(1,natoms+1):
		xnew=mol[i,'x']-x0
		ynew=mol[i,'y']-y0
		znew=mol[i,'z']-z0
		mol[i,'x']=xnew
		mol[i,'y']=ynew
		mol[i,'z']=znew


#------------------------------------------------------------
# function ztranslation
# translate ALL atoms in molecule
# ----------------------------------------------------------
def ztranslation(mol,distance):
	natoms=mol['natoms']
# translate all atoms	
	for i in range(1,natoms+1):
		znew=mol[i,'z']+distance
		mol[i,'z']=znew


#------------------------------------------------------------
# function ytranslation
# translate ALL atoms in molecule
# ----------------------------------------------------------
def ytranslation(mol,distance):
	natoms=mol['natoms']
# translate all atoms	
	for i in range(1,natoms+1):
		ynew=mol[i,'y']+distance
		mol[i,'y']=ynew


#------------------------------------------------------------
# function xtranslation
# translate ALL atoms in molecule
# ----------------------------------------------------------
def xtranslation(mol,distance):
	natoms=mol['natoms']
# translate all atoms	
	for i in range(1,natoms+1):
		xnew=mol[i,'x']+distance
		mol[i,'x']=xnew



#------------------------------------------------------------
# function bump_check_two_mols_flag1_under
# check for bumps between mo1l and mol2
# check only atoms where atom attribute under input value
# NB Want to exclude the two atoms forming the putative 
# bond between mol1 and mol2 
# These are atoms j1 in mol1 and j2 in mol2
# ----------------------------------------------------------
def bump_check_two_mols_flag1_under(value,mol1,mol2):
	natoms1=mol1['natoms']
	natoms2=mol2['natoms']
	dsqmax=10000
	bumpdist=sqrt(dsqmax)
	##print "n1 n2 = ",natoms1,natoms2
	for i in range(1,natoms1+1):
		for j in range(1,natoms2+1):
			flag1mol1=mol1[i,'flag1']
			flag1mol2=mol2[j,'flag1']
# Only check non deleted atoms
# that will be present in the combined molecule
# these have atom attribute 'flag1' less than value
			if (flag1mol1<value and flag1mol2<value):
				x1=mol1[i,'x']
				y1=mol1[i,'y']
				z1=mol1[i,'z']
				x2=mol2[j,'x']
				y2=mol2[j,'y']
				z2=mol2[j,'z']				 
				distsq=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
				if (distsq < dsqmax):
					dsqmax=distsq
					bumpdist=sqrt(distsq)
					##print "New shortest bump distance "
					##print "Distance = ",bumpdist
					##listatom(mol1,i)
					##listatom(mol2,j)
					##print "atoms numbers",i,a1,j,a2
					##print "-----------------------------------"
	return(bumpdist)


