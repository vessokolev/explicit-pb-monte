#!/usr/bin/env python

import numpy
import sqlite3
#
# Supply here the relative charge on each electrode - the number of electrons
# equivalent on each electrode with the respective sign:
#
left_electrode_charge=0.23
right_electrode_charge=-0.23
#
# Enter the number of grid nodes in x and y direction (nodes per Anstrom):
#
electrode_grid_density_x=0.2
electrode_grid_density_y=0.2
#
# The left electrode is at z=0. The right one is at the box size in z.
# On need to specify the core distance between the electrode and the ions
# so ions cannot approach the electrode bellow that distance (in Angstroms):
#
core_distance_at_electrode=3.0 # the gap between the electrode and the solution
#
# Core distance between the ions and the membrane in Angstroms.
#
core_distance_at_membrane=5.0 # the gap between the membrabe and the solution
#
# Set the PDB path and file name containing the input structure:
#
pdb_file="box_with_grid.pdb"
#
# Create version of the PDB file containing only the ions and save it in
# the file given bellow:
#
ion_pdb_file="ions.pdb"
#
# Set the path and filename of the file containing the produced Monte Carlo
# configuration:
#
inp_file="mc.inp"
#
# Supply the number of required successful proposals for Monte Carlo
#
accepted_proposals=30000000
#
# Number of box layers to get the long-range electrostatic interactions
#
number_of_box_layers=16
#
# Lennard-Jones cut-off parameter (in Angstroms):
#
LJ_cutoff=8.0


###############################################################
#             DO NOT MODIFY THE CODE BELLOW!!!!               #
###############################################################
#
# Create SQLite database and store in into the memory:
connection=sqlite3.connect(':memory:')
cursor=connection.cursor()
#
# Create a table for storing there the PDB atomic description
#
regexp='CREATE TABLE atoms(id INT PRIMARY KEY, serial INT, resName VARCHAR(3), \
        atomName VARCHAR(4), x REAL, y REAL, z REAL, sigma REAL, eps REAL, \
        spatial_g INT,charge REAL)'
cursor.execute(regexp)
#
# Create also the corresponding indexes to support the fast search:
#
regexp='CREATE INDEX atoms_idx ON atoms(resName,atomName,x,y,z)'
cursor.execute(regexp)
#
# van der Waals interaction parameters sigma and epsilon:
#
params=(\
["NA",3.3284,0.0115897],\
["CL",4.40104,0.4184])
#
# Relative charges of the atoms:
#
charge=(\
["NA",1.0],
["CL",-1.0])
#
# Reading the atoms from PDB file. Note that the description of each atom
# ther restart with the work ATOM or HETATM. Getting also the box dimensions
# from the header of the file.
#
counter=1

for line in open(pdb_file,"r"):
   check=line.split()
   if check[0]=='CRYST1':
      box=numpy.array([numpy.float(line[6:15]),\
                       numpy.float(line[15:24]),\
                       numpy.float(line[24:33])])
   if check[0]=='ATOM' or check[0]=='HETATM':
      name=line[0:6]
      serial=numpy.int(line[6:11])
      aname=line[12:16]
      resname=line[17:20]
      x=numpy.float(line[30:38])
      y=numpy.float(line[38:46])
      z=numpy.float(line[46:54])
      regexp='INSERT INTO atoms VALUES(?,?,?,?,?,?,?,?,?,?,?)'
      cursor.execute(regexp,(counter,serial,resname,aname,x,y,z,0,0,0,0,))
      counter+=1
#
# Create the ions only PBD file
#
f_obj=open(ion_pdb_file,"w")
for line in open(pdb_file,"r"):
   check=line.split()
   if check[0]=='ATOM' or check[0]=='HETATM':
      if check[3]=="NA" or check[3]=="CL":
         f_obj.write(line)
   else:
      f_obj.write(line)
f_obj.close()
#
# Get the both sides of the membrane in z
#
z_min=numpy.float(0.0)
z_max=numpy.float(0.0)

regexp='SELECT MIN(z) FROM atoms WHERE atomName=?'
cursor.execute(regexp,(" GRA",))

z_min=numpy.float(cursor.fetchall()[0][0])

regexp='SELECT MAX(z) FROM atoms WHERE atomName=?'
cursor.execute(regexp,(" GRA",))

z_max=numpy.float(cursor.fetchall()[0][0])

for i in params:
   regexp='UPDATE atoms SET sigma=? WHERE atomName=?'
   cursor.execute(regexp,(i[1],"  "+i[0]))
   regexp='UPDATE atoms SET eps=? WHERE atomName=?'
   cursor.execute(regexp,(i[2],"  "+i[0]))

for i in charge:
   regexp='UPDATE atoms SET charge=? WHERE atomName=?'
   cursor.execute(regexp,(i[1],"  "+i[0]))

regexp='UPDATE atoms SET spatial_g=1 WHERE z<?'
cursor.execute(regexp,(z_min,))

regexp='UPDATE atoms SET spatial_g=2 WHERE z>?'
cursor.execute(regexp,(z_max,))

#
# Open the output file for writing:
#
f_obj=open(inp_file,"w")
#
# Format the box sizes output:
#
line=""
line+="%9.3f" % box[0]
line+="%9.3f" % box[1]
line+="%9.3f" % box[2]
#
# Write down the box sizes:
#
f_obj.write(line+'\n')
#
# Write down the core distances:
#
line=""
line+="%12.8f" % core_distance_at_electrode
line+="%12.8f" % core_distance_at_membrane
f_obj.write(line+'\n')
#
# The same for the membrane charge:
#
line=""
line+="%12.6f" % left_electrode_charge
line+="%12.6f" % right_electrode_charge
f_obj.write(line+'\n')
#
# The electrode grid density (nodes per A)
#
line=""
line+="%8.3f" % electrode_grid_density_x
line+="%8.3f" % electrode_grid_density_y
f_obj.write(line+'\n')
#
# Do the same for the required number of accepted proposals:
#
line="%12d" % accepted_proposals
f_obj.write(line+'\n')
#
# Do the same for the required number of number of box images:
#
line="%5d" % number_of_box_layers
f_obj.write(line+'\n')
#
# Do the same for the Lennard-Jones cut-off parameter:
#
line="%12.8f" % LJ_cutoff
f_obj.write(line+'\n')
#
# Write down only the ions into separate PDB file to use it as a
# reference structure when visualize the trajectories.
#
regexp='SELECT * FROM atoms WHERE atomName="  NA" or atomName="  CL"'
cursor.execute(regexp)
result=cursor.fetchall()
#
# Start describing the atoms by supplying their total number:
#

line="%5d" % len(result)
f_obj.write(line+'\n')

counter=1

for i in result:
   line=""
   line+="%5d" % counter
   line+="%3s" % i[2]
   line+="%4s" % i[3]
   line+="%8.3f" % i[4]
   line+="%8.3f" % i[5]
   line+="%8.3f" % i[6]
   line+="%12.8f" % i[7]
   line+="%12.8f" % i[8]
   line+="%5d" % i[9]
   line+="%8.3f" % i[10]
   f_obj.write(line+'\n')
   counter+=1

f_obj.close()
