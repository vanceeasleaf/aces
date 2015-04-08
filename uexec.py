import sys
filename=sys.argv[1]
f=open(filename,"r")
fo=open("pbslist","w")
lmps=f.readlines()
for ll in lmps:
	