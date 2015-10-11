import numpy as np
from aces.tools import exists
class lineManager:
	def __init__(self,filename):
		self.f=open(filename)
		self._line=self.parse(self.f,filename)		
		self.nline=len(self._line)
	def parse(self,f,filename):
		npy=filename+'.npy'
		if exists(npy):return np.load(npy)
		print "scanning the file %s"%filename
		_line=[0]
		while f.readline():
			_line.append(f.tell())
		del _line[-1]
		np.save(npy,_line)
		return _line
		
	def getLine(self,i):
		self.moveto(i)
		return self.nextLine()
		
	def nextLine(self):
		return self.f.readline().strip()
	def moveto(self,i):
		f=self.f
		assert i < self.nline and i >=0
		f.seek(self._line[i])