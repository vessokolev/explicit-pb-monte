#!/usr/bin/env python

#
# This script is used to generate quasi-randomly distributed ions inside
# two neighboring boxes separated by a protein membrane, represented as
# a grid. There are two types of gaps that should be taken into account:
#
# 1. Next to the electrode planes (left one is at z=0 and the right one
#    is at the other end of the box in z). The role of the gap is to keep
#    away the ions from the box and protects this way the expressions for
#    computing of the Coulomb interactions from division by zero. The
#    length of the gap should be close to the adsorbtion layer thickness
#    in the Stern layer approach (sometimes addressed as Helmoltzh layer).
#    This way the simulation takes into account only the diffusion part of
#    the layer.
#
# 2. Next to the protein membrane. The ions should not penetrate the planes
#    of the membrane and should not come closer than ~ 5 Angstroms. 
# 
#
import numpy
import spread_atoms_across_box

###### DEFINE HERE THE PARAMETERS OF THE BOX ######

# This sizes (dimensions) of the simulation box:
box_x=40.0    # in x-direction (in Anstroms)
box_y=40.0    # in y-direction (in Anstroms)
box_z=170.0   # in z-direction (in Anstroms)

d=4.0 # grid spacing in Angstroms - the distance in x and y between the
      # neihboring nodes of the proteine membrane. It is used to generate
      # the positions of the dummy atoms of the protein membrane! It is NOT
      # involved with the ion-protein gap.

d_crit_wall=3.0 # The minimum distance between the wall (electrodes) and
                # the atoms (in Angstroms) in z-direction

d_gap_electrode=3.0 # The gap between the ions in the solution and the
                    # electode surface in Angstroms, in z-direction. The same
                    # gap is applied to the left and right electrodes. The size
                    # of the gap has the physical dimension of the Stern Layer.

d_gap_protein_grid=5.0 # The minimum distance between the outer layers of the
                       # protein grid and the ions in the solution (in A), in
                       # z-direction. An ion cannot get closer to the protein
                       # membrane than d_gap_protein_grid. (Number 5.0 is an
                       # approximation. There is no formula to compute it).

number_of_the_grid_plane_in_z=5 # How many planes (in z-direction) the grid
                                # have to have. Note that the protein membrane
                                # consists of perpendicular (in z) planes. Each
                                # plane contains on it certain amount of dummy
                                # atoms. Their surface density is 1/d/d.

# WARNING. If the protein grid membrane thickness is known in advance (in Angstroms)
# the number_of_the_grid_plane_in_z should be computed as:
#
# number_of_the_grid_plane_in_z=int(thickness/d)


# The path and name of the PDB file where the coordinates of the ions and grid
# atoms are stored:

output_file="box_with_grid.pdb"

C=0.5 # The concentration of the electrolite in the solution in M. It is used to
      # compute the number of Na and Cl ions.



###################################################
###### DO NOT EDIT THE CODE BELLOW THIS LINE ######
###################################################

repetitions_of_the_grid_plane_in_z=number_of_the_grid_plane_in_z

lower_x=0.0
lower_y=0.0
lower_z=0.0

# Keep the correct floating point format:
box_x=numpy.float(box_x)
box_y=numpy.float(box_y)
box_z=numpy.float(box_z)

#box_z-=2*numpy.float(d_gap_electrode)

lower_x=numpy.float(lower_x)
lower_y=numpy.float(lower_y)
lower_z=numpy.float(lower_z)

d=numpy.float(d)

# Angles between the unit vectors of the simulation box - 90 deg for prismatic box
ang_x=90.0
ang_y=90.0
ang_z=90.0

grid_node_x=int(box_x/d)+1
grid_node_y=int(box_y/d)+1

grid_nodes_2D_x=numpy.zeros(grid_node_x,dtype=numpy.float)
grid_nodes_2D_y=numpy.zeros(grid_node_y,dtype=numpy.float)
grid_nodes_2D_z=numpy.zeros(repetitions_of_the_grid_plane_in_z,dtype=numpy.float)


for p,q,r in zip([lower_x,lower_y],[grid_node_x,grid_node_y],[grid_nodes_2D_x,grid_nodes_2D_y]):
   r[0]=p
   for i in xrange(1,q):
      r[i]=r[i-1]+d

grid_nodes_2D_z[0]=lower_z
for i in xrange(1,repetitions_of_the_grid_plane_in_z):
   grid_nodes_2D_z[i]=grid_nodes_2D_z[i-1]+d

# Center the 2D nodes plane with respect to xy-plane
grid_nodes_2D_x+=(box_x-(grid_nodes_2D_x.max()-grid_nodes_2D_x.min()))/2.0
grid_nodes_2D_y+=(box_y-(grid_nodes_2D_y.max()-grid_nodes_2D_y.min()))/2.0
grid_nodes_2D_z+=(box_z-(grid_nodes_2D_z.max()-grid_nodes_2D_z.min()))/2.0


def write_atom_line(file_obj,atom_fmt,atom_vars):

   line=""
   for i,j in zip(atom_fmt,atom_vars):
      line+=i % j

   line+='\n'

   file_obj.write(line)

   return True


def write_cryst1_line(file_obj,cryst1_fmt,cryst1_vars):

   line=""

   for i,j in zip(cryst1_fmt,cryst1_vars):
      line+=i % j

   line+='\n'

   file_obj.write(line)

   return True


# The format of the CRYST1 PDB entry (this is the record that have to be placed in
# the first line of the PDB file). The CRYST1 entry allows to define the simulation
# box sizes.

lstring="P 21 21 21"
z_value=8

cryst1_fmt=["%6s","%9.3f","%9.3f","%9.3f","%7.2f","%7.2f","%7.2f","%1s","%11s","%4d"]
cryst1_vars=["CRYST1",box_x,box_y,box_z,ang_x,ang_y,ang_z,"",lstring,z_value]

# Open the output_file for writing:
file_obj=open(output_file,"w")
# Write the CRYST1 definitions in the output_file:
write_cryst1_line(file_obj,cryst1_fmt,cryst1_vars)
# Start writing the atoms in outpuf_file. Initialize the atomSerial number so the
# atom numbering starts from 1:
atomSerial=1
# Define the format vector for each ATOM line that is going to be written in the
# output_file:
atom_fmt=["%-6s","%5d","%1s","%4s","%1s","%3s","%1s","%1s","%4d","%1s","%3s","%8.3f","%8.3f","%8.3f","%6.3f","%6.3f","%10s","%2s","%2s"]
# Assign some atom name to the atoms of the grid. The name does not matter here
# but choose a distinguishable name:
atomName="GRA"
chainID=""
altLoc=""
resInsCode=""
occ=0.0
tempFact=0.0
atomCharge=""
elemSymb="G"
resName="GRD"
resNum=1



for i in xrange(grid_node_x):
   for j in xrange(grid_node_y):
      for q in xrange(repetitions_of_the_grid_plane_in_z):
         x=grid_nodes_2D_x[i]
         y=grid_nodes_2D_y[j]
         z=grid_nodes_2D_z[q]
         atom_vars=["ATOM",atomSerial,"",atomName,altLoc,resName,"",chainID,resNum,resInsCode,"",x,y,z,occ,tempFact,"",elemSymb,atomCharge]
         write_atom_line(file_obj,atom_fmt,atom_vars)
         atomSerial+=1


# Now lets guess the coordinates of the ions. The first task is to find their total
# number with respect to the concentration and the volume.

new_z=lower_z+grid_nodes_2D_z.min()-d_gap_protein_grid-d_gap_electrode

num_atoms=int(C*box_x*box_y*new_z*1.2046e-3)

if num_atoms % 2 <> 0:
   num_atoms+=1

cursor=spread_atoms_across_box.spread_atoms(box_x,box_y,new_z,num_atoms,0.5*new_z,0.8)

regexp='SELECT * FROM lattice'
cursor.execute(regexp)

counter=1

num_atoms/=2

ions_coords=cursor.fetchall()

shift_z=lower_z+d_crit_wall

for i in ions_coords:
   if counter<=num_atoms:
      atomName="NA"
      resName="NA"
      resNum=resNum+1
      x=i[0]
      y=i[1]
      z=i[2]+shift_z
      atom_vars=["ATOM",atomSerial,"",atomName,altLoc,resName,"",chainID,resNum,resInsCode,"",x,y,z,occ,tempFact,"",elemSymb,atomCharge]
      write_atom_line(file_obj,atom_fmt,atom_vars)
      atomSerial+=1
   else:
      atomName="CL"
      resName="CL"
      resNum=resNum+1
      x=i[0]
      y=i[1]
      z=i[2]+shift_z
      atom_vars=["ATOM",atomSerial,"",atomName,altLoc,resName,"",chainID,resNum,resInsCode,"",x,y,z,occ,tempFact,"",elemSymb,atomCharge]
      write_atom_line(file_obj,atom_fmt,atom_vars)
      atomSerial+=1
   counter+=1

del(cursor)

num_atoms*=2

cursor=spread_atoms_across_box.spread_atoms(box_x,box_y,new_z,num_atoms,0.5*new_z,2.0)

regexp='SELECT * FROM lattice'
cursor.execute(regexp)

counter=1

num_atoms/=2

ions_coords=cursor.fetchall()

shift_z=grid_nodes_2D_z.max()+d_gap_protein_grid

for i in ions_coords:
   if counter<=num_atoms:
      atomName="NA"
      resName="NA"
      resNum=resNum+1
      x=i[0]
      y=i[1]
      z=i[2]+shift_z
      atom_vars=["ATOM",atomSerial,"",atomName,altLoc,resName,"",chainID,resNum,resInsCode,"",x,y,z,occ,tempFact,"",elemSymb,atomCharge]
      write_atom_line(file_obj,atom_fmt,atom_vars)
      atomSerial+=1
   else:
      atomName="CL"
      resName="CL"
      resNum=resNum+1
      x=i[0]
      y=i[1]
      z=i[2]+shift_z
      atom_vars=["ATOM",atomSerial,"",atomName,altLoc,resName,"",chainID,resNum,resInsCode,"",x,y,z,occ,tempFact,"",elemSymb,atomCharge]
      write_atom_line(file_obj,atom_fmt,atom_vars)
      atomSerial+=1
   counter+=1


file_obj.close()

