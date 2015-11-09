## ## ## qdel.py v.0.1
## ## ## This program performs qdel for a range of values
## ## ## Created: 02/22/2013 - KDP

import os
import sys

assert len(sys.argv) == 3, "Must have a start and a stop argument"

start = sys.argv[1]
finish = sys.argv[2]

for i in range(start, finish):
	os.system('qdel '+ i)


