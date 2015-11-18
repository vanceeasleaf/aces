# encoding : utf8
import os,sys
import shutil
import pandas as pd
import subprocess as sub 
import shlex
import json
def loadjson(file):
    return json.loads(read(file))
def shell_exec(cmd):
	print "[Command]"+cmd
	sys.stdout.flush()
	c=os.popen(cmd).read()
	return c.strip()
def toString(m,sep=' '):
    return sep.join(map(str,m))
def passthru(cmd):
	#print os.popen(cmd).read()
	#sys.stdout.flush()
    #sub.call(shlex.split(cmd),stdout=sys.stdout)
	print "[Command]"+cmd
	sys.stdout.flush()
	sub.call(cmd,shell=True,stdout=sys.stdout)
def write(cmd,fileName,mode="w",sep=""):
	file=open(fileName,mode);
	file.write(str(cmd)+sep);
	file.close();
def debug(cmd):
    write(str(cmd),'debug.txt','a',sep="\n")
def sleep(u):
    os.sleep(u)
def exists(path):
    return os.path.exists(path)    
def read(fileName):
	file=open(fileName);
	s=file.read()
	file.close();
	return s
def to_txt(columns,data,filename):
	quants=pd.DataFrame(data,columns =columns)
	quants.to_csv(filename,sep='\t',index=None)
def pwrite(fp,s):
	print s,
	fp.write(s)
	sys.stdout.flush()	
def exit(info):
	print info
	sys.stdout.flush()
	sys.exit();
def mkcd(path):
    mkdir(path)
    cd(path)
def cd(path):
	os.chdir(path)
def mv(src,dest):
	shell_exec("mv %s %s"%(src,dest))
def cp(src,dest):
	shell_exec("cp %s %s -r "%(src,dest))
def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False
def ls(path='*'):
	import glob
	return glob.glob(path)
def pwd():
	return os.getcwd()
def dirname(path):
    return os.path.dirname(path).strip()
def parseyaml(filename):
    try:
        import yaml
    except ImportError:
        print "You need to install python-yaml."
        exit(1)
        
    try:
        from yaml import CLoader as Loader
        from yaml import CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper
    string = open(filename).read()
    data = yaml.load(string, Loader=Loader)
    return data
class File(file):
    """ An helper class for file reading  """

    def __init__(self, *args, **kwargs):
        super(File, self).__init__(*args, **kwargs)
        self.BLOCKSIZE = 4096

    def head(self, lines_2find=1):
        self.seek(0)                            #Rewind file
        return [super(File, self).next() for x in xrange(lines_2find)]

    def tail(self, lines_2find=1):  
        self.seek(0, 2)                         #Go to end of file
        bytes_in_file = self.tell()
        lines_found, total_bytes_scanned = 0, 0
        while (lines_2find + 1 > lines_found and
               bytes_in_file > total_bytes_scanned): 
            byte_block = min(
                self.BLOCKSIZE,
                bytes_in_file - total_bytes_scanned)
            self.seek( -(byte_block + total_bytes_scanned), 2)
            total_bytes_scanned += byte_block
            lines_found += self.read(self.BLOCKSIZE).count('\n')
        self.seek(-total_bytes_scanned, 2)
        line_list = list(self.readlines())
        return line_list[-lines_2find:]

    def backward(self):
        self.seek(0, 2)                         #Go to end of file
        blocksize = self.BLOCKSIZE
        last_row = ''
        while self.tell() != 0:
            try:
                self.seek(-blocksize, 1)
            except IOError:
                blocksize = self.tell()
                self.seek(-blocksize, 1)
            block = self.read(blocksize)
            self.seek(-blocksize, 1)
            rows = block.split('\n')
            rows[-1] = rows[-1] + last_row
            while rows:
                last_row = rows.pop(-1)
                if rows and last_row:
                    yield last_row
        yield last_row
