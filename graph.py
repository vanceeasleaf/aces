# encoding : utf8
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl	
def twinx(x,y1,y2,filename):
	fig,ax1=pl.subplots()
	a=ax1.plot(x[0],y1[0])
	ax1.set_xlabel(x[1])
	ax1.set_ylabel(y1[1])
	ax2=ax1.twinx()
	ax2.set_ylabel(y2[1])
	b=ax2.plot(x[0],y2[0],'r')
	ax1.set_frame_on(False)
	ax1.legend(a+b,[y1[1],y2[1]],loc='best').get_frame().set_alpha(0.0)
	pl.savefig(filename,bbox_inches='tight',transparent=True)

def plot(x,y,filename):
	series(x[1],y[1],[(x[0],y[0],y[1])],filename,legend=False)

def series(xlabel,ylabel,datas,filename,linewidth=1,legend=True):
	pl.figure()
	pl.xlabel(xlabel)
	pl.ylabel(ylabel)
	for serie in datas:
		pl.plot(serie[0],serie[1],label=serie[2],linewidth=linewidth)
	if legend:
		pl.legend(loc='best').get_frame().set_alpha(0.0)
	pl.savefig(filename,bbox_inches='tight',transparent=True) 