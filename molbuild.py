#!/usr/bin/env python3


import sys
import os
# access functions in mol_utils.py file
from mol_utils import *
from math import *

#test exec and system calls


#------------------------------------------------------------
# function check_file_exists
# ERROR and stop program if named file does not exist
# ------------------------------------------------------------
def check_file_exists(fname):
	if (not os.path.isfile(fname)):
		print("ERROR file not found ",fname)
		print("In script molbuild.py")
		print("Program must exit ")
		sys.exit(1)

#------------------------------------------------------------
# function command_use
# Processing for the 'use' command
# example 1:
# use bdglc.pdb as mol1
# which results in an OS command to copy bdglc.pdb to temp_mol1.pdb
# example 2: 
# use frag_4adglc as mol2 resname ADG resumber 2 
# copies input file as above but runs residue renaming/renumbering script
# NB resname and resnumber features not added yet !!!!
# ------------------------------------------------------------
def command_use(command):
	if (command.find("use")<0):
		print("String 'use' not found for this command",command)
		print("Program must exit ")
		sys.exit(1)		
	# the 'use' command was found
	## print "Use command found "
	# create a list holding each word
	words=command.split()
	#print "words =",words[0],words[1],words[2],words[3]
	#print "words 0=",words[0]
	#print "words 1=",words[1]
	#print "words 2=",words[2]
	#print "words 3=",words[3]
	if (not(words[0]=="use")):
		print("First word of this command should be 'use', command=",command)
		print("Program must exit ")
		sys.exit(1)		
	if (not(words[2]=="as")):
		print("Third word of this command should be 'as', command=",command)
		print("Program must exit ")
		sys.exit(1)

	# words[1] must be a valid file or path/file name
	# function will exit program if an error found
	check_file_exists(words[1])
	# The command - 'use bdglc.pdb as mol1'
	# should result in an OS command 'cp bdglc.pdb temp_mol1.pdb'
	# hence words[1] should contain a file name such as bdglc.pdb
	#if (os.path.isfile(words[1])):
	pdbfile=words[1]
	##print "using file ",pdbfile
	# fourth word words[3] used to construct working file name
	newpdbfile="temp_"+words[3]+".pdb"
	## print "new file=",newpdbfile
	oscmd="cp "+pdbfile+" "+newpdbfile		
	## print "oscmd=",oscmd						
	os.system(oscmd)

#	now we check if 'resname' or 'resnumber' commands found
#	as part of the 'use' command
#	example command 'use bdglc.pdb as mol1 resname BDG resnumber 2'
#	the resname and resnumber command will operate on file (eg) temp_mol1.pdb
#
	resname_found=words.count('resname')
	resnumber_found=words.count('resnumber')
	
	if (resname_found):
		resname_index=words.index('resname')
		new_resname=words[1+resname_index]
		#print "Set new residue name to ",new_resname
		resname_str="resname "+newpdbfile+" "+new_resname
		#print "residue name command to execute =",resname_str
		command_resname(resname_str)

	if (resnumber_found):
		resnumber_index=words.index('resnumber')
		new_resnumber=words[1+resnumber_index]
		#print "Set new residue number to ",new_resnumber
		resnumber_str="resnumber "+newpdbfile+" "+new_resnumber
		#print "residue number command to execute =",resnumber_str
		command_resnumber(resnumber_str)




#------------------------------------------------------------
# function command_resname
# Processing for the 'resname' command
# example 1:
# resname mol1 GLC
# which leads to OS command
# pdb_set_resname.py temp_mol1.pdb GLC
# this function creates file temp_temp.pdb and copies to temp_mol1.pdb
# File temp_mol1.pdb must have been created
# by a previous part of 'use' command
# NB Now called from command_use
# ------------------------------------------------------------
def command_resname(command):
	if (command.find("resname")<0):
		print("String 'resname' not found for this command",command)
		print("Program must exit ")
		sys.exit(1)	
	# residue name command e.g.
	# pdb_set_resname.py mol1 GLC
	# where mol1 must have been created as file temp_mol1.pdb
	# by using a 'use' command as above
	words=command.split()
	if (words[0]=="resname"):
	## print "words 0 1 2 = ",words[0],words[1],words[2]
		##molfilename="temp_"+words[1]+".pdb"
		molfilename=words[1]
		check_file_exists(molfilename)
		resname=words[2]
		resname2=resname.upper()
		oscmd1="./pdb_set_resname.py"+" "+molfilename+" "+resname2+" > temp_temp.pdb"
		oscmd2="cp temp_temp.pdb "+molfilename
		#print "RESNAME OSCMD 1 = ",oscmd1
		#print "RESNAME OSCMD 2 = ",oscmd2
		os.system(oscmd1)
		os.system(oscmd2)

#------------------------------------------------------------
# function command_save
# Processing for the 'save' command
# example 1:
# save keep_mol1.pdb
# which leads to OS command
# to copy temp_current.pdb to keep_mol1.pdb
# NB - ONLY use after 'join' commands
# ------------------------------------------------------------
def command_save(command):
	if (command.find("save")<0):
		print("String 'save' not found for this command",command)
		print("Program must exit ")
		sys.exit(1)
	words=command.split()
	if (words[0]=="save"):
		molfilename=words[1]
		oscmd1="cp temp_current.pdb "+molfilename
		#print "SAVE OSCMD  = ",oscmd1
		os.system(oscmd1)

#------------------------------------------------------------
# function command_resnumber
# Processing for the 'resnumber' command
# example 1:
# resnumber mol1 1
# which leads to OS command
# pdb_set_resnumber.py temp_mol1.pdb 1
# where mol1 must have been created as file temp_mol1.pdb
# by a previous 'use' part of command
# Now called from command_use function
# ------------------------------------------------------------
def command_resnumber(command):
	if (command.find("resnumber")<0):
		print("String 'resnumber' not found for this command",command)
		print("Program must exit ")
		sys.exit(1)

	words=command.split()
	if (words[0]=="resnumber"):
		##print "words 0 1 2 = ",words[0],words[1],words[2]
		##molfilename="temp_"+words[1]+".pdb"
		molfilename=words[1]		
		check_file_exists(molfilename)
		resnum=words[2]
		resnum2=resnum.upper()
		oscmd1="./pdb_set_resnumber.py"+" "+molfilename+" "+resnum2+" > temp_temp.pdb"
		oscmd2="cp temp_temp.pdb "+molfilename
		#print "RESNUMBER OSCMD 1 = ",oscmd1
		#print "RESNUMBER OSCMD 2 = ",oscmd2
		os.system(oscmd1)
		os.system(oscmd2)


#------------------------------------------------------------
# function command_join
#
# join command
# examples:
# join mol1 1.O1 1.C1 fragment mol2
# join current 2.O2 2.C2 fragment mol3
#
# 
# ------------------------------------------------------------
def command_join(command):
	if (command.find("join")>=0):
		words=command.split()
		if (words[0]=="join"):
			#print "words 0 1 2 3 4 5= ",words[0],words[1],words[2],words[3],words[4],words[5]
			mol1filename="temp_"+words[1]+".pdb"
			mol2filename="temp_"+words[5]+".pdb"
			check_file_exists(mol1filename)
			check_file_exists(mol2filename)

#
# files temp_mol1.pdb and temp_mol2.pdb do exist - proceed
#
			delatom=words[2]
			joinatom=words[3]
#
# could check words[4] is 'fragment' but not doing so yet
#
			# if command of form
			# join current 3.O1 3.C1 fragment mol4
			# then need special action so that file temp_current is NOT overwritten by redirection
			# hence make a copy to temp_temp.pdb and use this as input molecule file
			# output then redirected to temp_current.pdb
			if (words[1]=="current"):
			# need file temp_current.pdb for output so copy to temp_temp.pdb
			# and adjust mol1filename - not mol1, mol2 etc.				
				oscmd2="cp temp_current.pdb temp_temp.pdb"
				#print "JOIN OSCMD =",oscmd2
				os.system(oscmd2)
				mol1filename="temp_temp.pdb"
			keepfile="temp_current.pdb"
			keepstring=" > "+keepfile
			oscmd1="./pdb_join_fragment.py"+" "+mol1filename+" "+delatom+" "+joinatom+" "+mol2filename+keepstring
			#print "JOIN OSCMD  ",oscmd1
			os.system(oscmd1)	



#------------------------------------------------------------
# function process_command
# process each line as a single command
# MAYBE add 'save' command later
# ------------------------------------------------------------
def process_command(line2):
	## print "Command =",line2
##
## the 'use' command
## example:
## use bdglc.pdb as mol1
##
	if (line2.find("use ")>=0):
		command_use(line2)
##
## join command
## examples:
## join mol1 1.O1 1.C1 fragment mol2
## join current 2.O2 2.C2 mol3
##

	if (line2.find("join ")>=0):
		command_join(line2)	
##
## Added 'save' command
## e.g. save hep4mer.pdb 
## copies temp_current.pdb to hep4mer.pdb
##
	if (line2.find("save ")>=0):
		command_save(line2)


#------------------------------------------------------------
# function read_commands
# read molecule building commands from file and execute them
# ------------------------------------------------------------
def read_commands(filename):
	print('Commands filename = ',filename)
	commandsfile=open(filename,'r')
	for line in commandsfile.readlines():
		print(line)
		line.lstrip
		##process_command(line)
# if line starts with a hash character then it is a comment and is ignored
		if (line[0]=="#"):
			comment_line=1
# could print comment to STDERR
		else:
			process_command(line)




# ------------------start of main program----------------------

# command line arguments
# programname=sys.argv[0]
#print 'Number of arguments = ',len(sys.argv)
if len(sys.argv)!=2:
	print('Usage molbuild.py instructions.txt')
	print('e.g. ./molbuild.py build_maltotriose.txt')
	sys.exit(1)
	
filename=sys.argv[1]

##print "filename =",filename
##command="cat "+filename
##print "command =",command
##os.system(command)
read_commands(filename)
print("Done")






