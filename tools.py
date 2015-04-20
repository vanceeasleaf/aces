# encoding : utf8
import os,sys
import shutil
def shell_exec(cmd):
	c=os.popen(cmd).read()
	return c.strip()
def passthru(cmd):
	print os.popen(cmd).read()
	sys.stdout.flush()	
def write(cmd,fileName):
	file=open(fileName,"w");
	file.write(cmd+"\n");
	file.close();
def exit(info):
	print info
	sys.stdout.flush()
	sys.exit();
def cd(path):
	os.chdir(path)
def mv(src,dest):
	shutil.move(src,dest)
def cp(src,dest):
	shutil.copyfile(src,dest)
def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False