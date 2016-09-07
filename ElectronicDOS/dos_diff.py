#!/usr/bin/env python

# dos_diff.py v1.0 9-18-2012 Jeff Doak jeff.w.doak@gmail.com

from electronicdos import *
import sys

dos1dir = str(sys.argv[1])
dos2dir = str(sys.argv[2])
dos1 = ElectronicDOS(dos1dir+"/DOSCAR",dos1dir+"/OUTCAR",dos1dir+"/POSCAR")
dos2 = ElectronicDOS(dos2dir+"/DOSCAR",dos2dir+"/OUTCAR",dos2dir+"/POSCAR")
#dos1.shift_energy(-1.*dos1.e_fermi)
if len(sys.argv) > 3:
    dos1.shift_energy(float(sys.argv[3]))
#else:
#    dos2.shift_energy(-1.*dos2.e_fermi)
spline1 = dos1.dos_spline()
spline2 = dos2.dos_spline()
difference = dos1.dos_difference(spline1,spline2)
for i in range(len(dos2.energy)):
    print dos2.energy[i],spline2(dos2.energy[i]),spline1(dos2.energy[i]),difference[i]
sys.exit()
