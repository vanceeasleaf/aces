class lineManager:
	def __init__(self,filename):
		self.f=open(filename)
		self._line=[0]
		self.parse()		
	def parse(self):
		f=self.f
		while f.readline():
			self._line.append(f.tell())
		del self._line[-1]
		self.nline=len(self._line)
		
	def getLine(self,i):
		self.moveto(i)
		return self.nextLine()
		
	def nextLine(self):
		return self.f.readline().strip()
	def moveto(self,i):
		f=self.f
		assert i < self.nline and i >=0
		f.seek(self._line[i])