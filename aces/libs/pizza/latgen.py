# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# latgen tool

oneline = "Convert LAMMPS snapshots to latgen format"

docstr = """
x = latgen(d)		d = object containing atom coords (dump, data)

x.one()                 write first snapshots to latgen
x.single(N)             write snapshot for timestep N to latgen 
x.single(N,"file")      write snapshot for timestep N to file.vasp
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   data = data file to read from

# Imports and external programs

import sys
import numpy as np
# Class definition

class latgen:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data
   
  # --------------------------------------------------------------------

  def one(self):
    self.single(0)

      
  
  
  # --------------------------------------------------------------------

  def single(self,time,*args):
    if len(args) == 0: file = "cell.in"
    #elif args[0][-5:] == ".vasp": file = args[0]
    else: file = args[0]
    #self.data.scale()
    which = self.data.findtime(time)
    time,box,atoms,bonds,tris,lines = self.data.viz(which)
    f = open(file,"w")
    #print >>f,self.data.title,
    print >>f,1.0 #lattice constant
    xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz=box[0],box[1],box[2],box[3],box[4],box[5],box[6],box[7],box[8]
    lx=xhi-xlo
    ly=yhi-ylo
    lz=zhi-zlo
    print >>f,"%f\t0.0\t0.0" % (lx)
    print >>f,"%f\t%f\t0.0" % (xy,ly)
    print >>f,"%f\t%f\t%f" % (xz,yz,lz)
    #print >>f,"Cartesian"
    #print >>f,len(atoms)
    atoms=self.scale(box,atoms)
    typesatom=[atom[1] for atom in atoms ]
    types=list(set(typesatom))
    ntype=len(types)
    atomOfType={}
    for type in types:
        atomOfType[type]=[atom for atom in atoms if atom[1]==type]
    for type in types:
        print >>f,"%d\t" % (len(atomOfType[type])),
    print >>f,""
    for type in types:
	  for atom in atomOfType[type]:
		  itype = int(atom[1])
		  print >>f,atom[2],atom[3],atom[4]

    f.close()
    
  def scale(self,box,atoms):
    xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz=box[0],box[1],box[2],box[3],box[4],box[5],box[6],box[7],box[8]
    atoms=np.array(atoms)
    if 0 and xy == 0.0 and xz == 0.0 and yz == 0.0:
      xprdinv = 1.0 / (xhi - xlo)
      yprdinv = 1.0 / (yhi - ylo)
      zprdinv = 1.0 / (zhi - zlo)

      atoms[:,2] = (atoms[:,2] - xlo) * xprdinv
      atoms[:,3] = (atoms[:,3] - ylo) * yprdinv
      atoms[:,4] = (atoms[:,4] - zlo) * zprdinv
    else:
      h0 = xhi - xlo
      h1 = yhi - ylo
      h2 = zhi - zlo
      h3 = yz
      h4 = xz
      h5 = xy
      h0inv = 1.0 / h0
      h1inv = 1.0 / h1
      h2inv = 1.0 / h2
      h3inv = yz / (h1*h2)
      h4inv = (h3*h5 - h1*h4) / (h0*h1*h2)
      h5inv = xy / (h0*h1)
      atoms[:,2] = (atoms[:,2] - xlo)*h0inv + \
          (atoms[:,3] - ylo)*h5inv + \
          (atoms[:,4] - zlo)*h4inv
      atoms[:,3] = (atoms[:,3] - ylo)*h1inv + \
          (atoms[:,4] - zlo)*h3inv
      atoms[:,4] = (atoms[:,4] - zlo)*h2inv
    return atoms
