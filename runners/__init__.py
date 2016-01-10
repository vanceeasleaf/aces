import sys
from aces.tools import shell_exec
class Runner:
	def __init__(self,m):
		self.m=m
	def generate(self):
		pass
	def run(self):
		__console__=sys.stdout
		
		f=open('input', 'w',0)
		sys.stdout=f
		self.generate()
		sys.stdout=__console__
		f.close()
		cmd=self.runcmd()
		shell_exec(cmd)
		self.post()
	def runcmd(self):
		return ""
	def post(self):
		pass