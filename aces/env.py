from os.path import *
from aces.tools import pwd,exists
import sys
SRCHOME=dirname(realpath(__file__))
def checkParent(dir,n=0):
	if n==5:raise Exception('error when find sub.py')
	if not exists(dir+'/sub.py'):
		dir=realpath(dir+'/..')
		return checkParent(dir,n+1)		
	else:
		return dir
PROJHOME=checkParent(pwd())
PROJNAME=basename(PROJHOME)