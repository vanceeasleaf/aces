import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl


def drawTempAve():
	f=open('tempAve.txt')
	f.readline()
	x=[];t=[];j=[]
	for line in f:
		id,Coord,Count,Temp,Jx=line.strip().split()
		x.append(Coord);
		t.append(Temp);
		j.append(Jx);
	f.close()
	fig,ax1=pl.subplots()
	#x=x[3:-3];t=t[3:-3];j=j[3:-3]
	a=ax1.plot(x,t)
	ax1.set_ylabel('temperature(K)')
	ax1.set_xlabel('x(Augstrom)')
	ax2=ax1.twinx()
	ax2.set_ylabel('heat flux')
	b=ax2.plot(x,j,'r')
	ax1.set_frame_on(False)
	ax1.legend(a+b,["temperature","heat flux"],loc='best').get_frame().set_alpha(0.0)
	pl.savefig('tempAve.png',bbox_inches='tight',transparent=True)
 
if __name__=='__main__':
	drawTempAve()