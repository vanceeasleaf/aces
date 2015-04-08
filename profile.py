#encoding : utf8
class profile:
	def getTempProfile(self,begin,upP,deta,S,tcfactor,zfactor):
		pass
		
if __name__=='__main__':
	import sys
	begin,upP,deta,S,tcfactor,zfactor=sys.argv[1:]
	profile().getTempProfile(int(begin),int(upP),float(deta),float(S),float(tcfactor),float(zfactor))