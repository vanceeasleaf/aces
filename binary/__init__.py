from aces.tools import passthru,cd,shell_exec
from os.path import *
def home(name):
	return dirname(realpath(__file__))+"/"+name
def run(name):
	path=home(name)
	if not exists(path):
		workdir=shell_exec('pwd')
		cd(dirname(path))
		passthru("make")
		cd(workdir)
	passthru(path)
def pr():
	run('pr')
if __name__=='__main__':
	pr()