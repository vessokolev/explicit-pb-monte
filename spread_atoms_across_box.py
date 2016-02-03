#!/usr/bin/env python

import numpy
import sqlite3

def spread_atoms(size_x,size_y,size_z,num_atoms,core_dist_init,core_dist_final):
   #
   # In order to use NumPy functions properly the variables should be
   # properly formated to fit the types of the arguments.
   #
   # size_x and size_y are the sizes of the planar area where the
   # quasi-random points should be placed.
   #
   # core_distance_sqr stands for the critical minimum distances between
   # the points on the plane.
   #
   size_x=numpy.float(size_x)
   size_y=numpy.float(size_y)
   size_z=numpy.float(size_z)
   core_distance_sqr=numpy.float(core_dist_init*core_dist_init)
   core_distance_sqr_low=numpy.float(core_dist_final*core_dist_final)

   numRandVectors=100000
   #
   # To implement Periodic Boundary Conditions (PBC) the halfs of each
   # size of the rectangular should be computed and saved as parameter for
   # using in the PBC selection routine. In principle it is possible to get
   # these values on the fly (on request) every time the algorithm needs
   # them but since they are supposed to be constant here there is no such
   # need. Also defining them as constants speeds up the execution.
   #
   size_x_half=0.5*size_x
   size_y_half=0.5*size_y
   size_z_half=0.5*size_z
   #
   # The SQLite database is created and maintained in the memory. This is
   # the fastest and most efficient way of manipulating small SQLIte
   # databases.
   #
   connection=sqlite3.connect(":memory:")
   cursor=connection.cursor()
   #
   # The table named "lattice" is a storage object in SQLite database
   # used for storing the x and y-coordinates of each accepted vector.
   # This table should also store the critical minimum distance between
   # two nodes, used for guessing the accepted vector, in a separate column.
   #
   regexp="CREATE TABLE lattice(x REAL, y REAL, z REAL, d REAL)"
   cursor.execute(regexp)
   #
   # To increase the speed of the seclection an index need to be created.
   # N1ote that the critical minimum distance paramter, d, should not be
   # indexed because its column is not involved in any selection done
   # by this code.
   #
   regexp="CREATE INDEX lattice_indx ON lattice(x,y,z)"
   cursor.execute(regexp)
   #
   # Define four 3D vectors which will be used to hold the reflections
   # of the probbing vector with respect to the PBC.
   #
   reflections=numpy.zeros(24,dtype=numpy.float).reshape(8,3)

   numHits=0
   insertedAtoms=0
   #
   # Trying to use the same d for testing all vectors can lead to
   # inifinity long simulation. But if the target is to generate a
   # quasirandomly distributed points it is possible to increment
   # the critical minimum distance, d. That process should takes place
   # when the number successively rejected random vector reaches
   # the value of max_tries.
   #
   max_tries=10000
   #
   # When the number of rejections becomes equal to max_tries the
   # critical minimum distance, d, should be incremented by
   # multiplying it by decrement:
   #
   # d=d*decrement
   #
   decrement=0.99
   counter=0
   #
   # The variable flag controls the first of the nested 'while' loops.
   # It have to be set True to start the loop. During the loop its
   # values might be set False if some of the criterions for canceling
   # the loop process is met.
   #
   flag=True

   while flag:

      randVectors=numpy.zeros(3*numRandVectors,dtype=numpy.float).reshape(numRandVectors,3)
      randVectors[:,0]=numpy.random.rand(numRandVectors)*size_x
      randVectors[:,1]=numpy.random.rand(numRandVectors)*size_y
      randVectors[:,2]=numpy.random.rand(numRandVectors)*size_z
      i=0
      randVectorsHighestIndex=numRandVectors-1

      while insertedAtoms<num_atoms or core_distance_sqr_low>core_distance_sqr:

            if i==randVectorsHighestIndex:
               randVectors[:,0]=numpy.random.rand(numRandVectors)*size_x
               randVectors[:,1]=numpy.random.rand(numRandVectors)*size_y
               randVectors[:,2]=numpy.random.rand(numRandVectors)*size_z
               i=0

         
            flag_x=randVectors[i,0]<=size_x_half
            flag_y=randVectors[i,1]<=size_y_half
            flag_z=randVectors[i,2]<=size_z_half
            #
            # Computing the PBC reflections of the probing vector:
            #
            if flag_x and flag_y:
               reflections[0]=[randVectors[i,0],randVectors[i,1]+size_y,randVectors[i,2]]
               reflections[1]=[randVectors[i,0]+size_x,randVectors[i,1]+size_y,randVectors[i,2]]
               reflections[2]=[randVectors[i,0]+size_x,randVectors[i,1],randVectors[i,2]]

            if flag_x and not flag_y:
               reflections[0]=[randVectors[i,0]+size_x,randVectors[i,1],randVectors[i,2]]
               reflections[1]=[randVectors[i,0]+size_x,randVectors[i,1]-size_y,randVectors[i,2]]
               reflections[2]=[randVectors[i,0],randVectors[i,1]-size_y,randVectors[i,2]]

            if not flag_x and not flag_y:
               reflections[0]=[randVectors[i,0],randVectors[i,1]-size_y,randVectors[i,2]]
               reflections[1]=[randVectors[i,0]-size_x,randVectors[i,1]-size_y,randVectors[i,2]]
               reflections[2]=[randVectors[i,0]-size_x,randVectors[i,1],randVectors[i,2]]

            if not flag_x and flag_y:
               reflections[0]=[randVectors[i,0]-size_x,randVectors[i,1],randVectors[i,2]]
               reflections[1]=[randVectors[i,0]-size_x,randVectors[i,1]+size_y,randVectors[i,2]]
               reflections[2]=[randVectors[i,0],randVectors[i,1]+size_y,randVectors[i,2]]

            reflections[3]=randVectors[i,:]

            if flag_z:
               for p,q in zip(xrange(4,8),xrange(4)):
                  reflections[p]=reflections[q]
                  reflections[p][2]+=size_z
            else:
               for p,q in zip(xrange(4,8),xrange(4)):
                  reflections[p]=reflections[q]
                  reflections[p][2]-=size_z

            for j in xrange(4):
               regexp="SELECT COUNT(*) FROM lattice WHERE ((?-x)*(?-x)+(?-y)*(?-y)+(?-z)*(?-z))<=?"
               cursor.execute(regexp,(reflections[j,0],reflections[j,0],reflections[j,1],reflections[j,1],reflections[j,2],reflections[j,2],core_distance_sqr,))
               numHits+=cursor.fetchall()[0][0]

            counter+=1

            if numHits==0:
               #
               # If no collision is found the probbing vector should be accepted and stored
               # in the SQLite database. Note that each vector added to the dabase gets its
               # x,y-coordinates indexed to support the search efficiency. The values of
               # critical minimum distance, 
               #
               regexp="INSERT INTO lattice VALUES(?,?,?,?)"
               cursor.execute(regexp,(randVectors[i,0],randVectors[i,1],randVectors[i,2],core_distance_sqr,))
               insertedAtoms+=1
#               print randVectors[i,0],randVectors[i,1],randVectors[i,2],core_distance_sqr
               print "New random vector has been found ..."
               counter=0

            if counter==max_tries:
               #
               #
               #
               core_distance_sqr*=decrement**2.0
               counter=0
            #
            # Make numHits zero to prepare it to enter the loop again (if needed).
            #
            numHits=0
   
            i+=1
   
      flag=False

#   connection.commit()
#   connection.close()

   return cursor



