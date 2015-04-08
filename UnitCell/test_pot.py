
import numpy as np
from unitcell import *
import re

#for i in range(len(lines)):
#    if lines[i].startswith("  (the norm of the test charge is"):
#        j = 1
#        k = 0
#        while k < poscar.num_atoms:
#            line = lines[i+j].split()
#            while len(line) > 0:
#                k = int(line.pop(0))
#                charged_pots[k-1] = float(line.pop(0))
#            j += 1
#        break



def get_el_pots(outcar_name,n_atoms):
    outcar = open(outcar_name,"r")
    lines = outcar.readlines()
    outcar.close()
    charged_pots = np.zeros(n_atoms)
    for i in range(len(lines)):
        if lines[i].startswith("  (the norm of the test charge is"):
            j = 1; k = 0;
            while k < n_atoms:
                line = lines[i+j].split()
                print "Original Line:",line
                while len(line) > 0:
                    print line
                    print
                    temp = line.pop(0)
                    print "Line is now:",line,len(line)
                    print "Trying to grab value from string:",temp
                    try:
                        k = int(temp)
                    except ValueError:
                        el_str = str(k+1)+r"([-][0-9][0-9]*[.][0-9][0-9]*)"
                        el_reg = re.compile(el_str)
                        k += 1
                        print el_str,temp
                        charged_pots[k-1] = float(el_reg.match(temp).group(1))
                        print k,charged_pots[k-1]
                    else:
                        charged_pots[k-1] = float(line.pop(0))
                        print "Found value:",k-1,charged_pots[k-1]
                j += 1
            break
    return charged_pots



chargedout="OUTCAR"
chargedpos="POSCAR"
poscar = UnitCell(chargedpos)
n_atoms = poscar.num_atoms
charged_pots = get_el_pots(chargedout,n_atoms)
print charged_pots
