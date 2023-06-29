# pdb_climb
PDB climb - command line interface molecule builder for molecular structures in protein databank format.

The code PDBclimb code constitutes a very basic command line interface molecule builder for molecules
in protein data bank (PDB) format.


It is intended as a set of scripts to be used as 'molecular scissors'. 
This code was created to meet the following basic requirements
	
	 * Snipping and joining' 3D molecular graphs - for molecules stored PDB format files.
	 * Specify the atoms involved by atom name and residue number - not a numeric atom index
	 * Allow very specific control over residue names and numbers as part of the building process

Note that this process occurs by 'snipping' and 'joining' single bonds, not double or triple bonds.
Note that the snipping and joining single bond should not be a ring bond. 

As the name implies it operates only on protein data bank (PDB) format files. The code is written 
in Python (was Python2 but now Python3) and should run easily on any machine with Python3 installed.


The scripts were created for the purpose of building PDB format files of small oligosaccharide 
sugar molecules, and many of these examples reflect this. However the codes are completely general 
and can be used for any type of molecule in PDB format. The code is relatively simple in that
it focuses on the so called HETATM and ATOM records of PDB files. It does not seek to manipulate 
and deal with the CONECT records specifying bonds between atoms. 


Please see the Markdown files in the documentation folder for an introduction to the code and for 
a description of the individual python scripts. 

Having downloaded this repository you should have the following folders. 

* documentation - This contains documentation in Markdown format and images used by these files.
* pdb - A folder of potentially usedful PDB format files. 
* frag - A folder of 'fragment' PDB files. These use the normal PDB format but follow the conventions 
		described in the introductory documentation.
* build - A series of molecule build files as used by the molbuild.py script. 

Hope this is useful.

