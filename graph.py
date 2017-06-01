# encoding : utf8
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.width'] = 1.5
matplotlib.rcParams['ytick.major.width'] = 1.5
matplotlib.rcParams['patch.linewidth']=0
matplotlib.rcParams['legend.markerscale']=.7
matplotlib.rcParams['mathtext.default']='regular'
#matplotlib.rcParams['legend.frameon'] = 'False'
import numpy as np
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
	pl.close()
def plot(x,y,filename,**args):
	args['legend']=False
	series(x[1],y[1],[(x[0],y[0],y[1])],filename,**args)
	
def series(xlabel,ylabel,datas,filename,linewidth=1,legend=True,grid=False,xmax=None,xmin=None,scatter=False,logx=False,logy=False):
	pl.figure()
	pl.xlabel(xlabel)
	pl.ylabel(ylabel)
	plot=pl.plot
	if logx:
		plot=pl.semilogx
	if logy:
		plot=pl.semilogy
	if logx and logy:
		plot=pl.loglog
	if scatter:
		marker='.'
		linestyle='.'
	else: 
		marker=None
		linestyle='-'
	if len(datas)==1:
		serie=datas[0]
		plot(serie[0],serie[1],label=serie[2],linewidth=linewidth,color='r',marker=marker,linestyle=linestyle)
		if xmax is None:xmax=max(serie[0])
		if xmin is None:xmin=min(serie[0])
		pl.xlim([xmin,xmax])
	else:
		min0=100000;max0=-100000
		for serie in datas:
			plot(serie[0],serie[1],label=serie[2],linewidth=linewidth,marker=marker,linestyle=linestyle)
			min0=min(min(serie[0]),min0)
			max0=max(max(serie[0]),max0)
		if xmax is None:xmax=max0
		if xmin is None:xmin=min0
		pl.xlim([xmin,xmax])
	if legend:
		pl.legend(loc='best').get_frame().set_alpha(0.0)
	if grid:
		pl.grid(True)
	
	pl.savefig(filename,bbox_inches='tight',transparent=True) 
	pl.close()
def wfig(filename,f,legend=False):
	pl.figure()
	f()
	
	if legend:
		pl.legend(loc='best').get_frame().set_alpha(0.0)
	pl.savefig(filename,bbox_inches='tight',transparent=True) 
	pl.close()
def setLegend(pl,ncol=1):
	frame=pl.legend(loc='best',scatterpoints=1,numpoints=1,ncol=ncol,fontsize=12).get_frame()
	frame.set_linewidth(1.5)
class fig:
	def __init__(self,filename,legend=False,ncol=3,figsize=None):
		self.filename=filename
		self.legend=legend
		self.ncol=ncol
		self.figsize=figsize
	def __enter__(self):
		pl.figure(figsize=self.figsize)
	def __exit__(self, type,value, trace):
		if self.legend:
			setLegend(pl,ncol=self.ncol)
		pl.savefig(self.filename,bbox_inches='tight',transparent=True) 
		pl.close() 
def surf(X, Y, Z,filename):
	from mpl_toolkits.mplot3d import Axes3D
	fig = pl.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X, Y, Z)
	pl.savefig(filename,bbox_inches='tight',transparent=True) 
	pl.close()
def imshow(c,filename,extent=[0,1,0,1]):
	pl.figure()
	pl.imshow(np.array(c).astype('float'),extent=extent,origin='lower')
	pl.savefig(filename,bbox_inches='tight',transparent=True) 
	pl.close()
def scatter(x,y,c,xlabel,ylabel,filename,marker_size=2.0,marker='s'):
	pl.figure()
	pl.xlabel(xlabel)
	pl.ylabel(ylabel)
	pl.scatter(x,y,edgecolors='none',s=marker_size,c=c,marker=marker)
	pl.savefig(filename,bbox_inches='tight',transparent=True) 
	pl.close()
def scatter3d(x,y,z,filename):
	from mpl_toolkits.mplot3d import Axes3D,proj3d
	fig = pl.figure()
	ax = fig.add_subplot(111, projection='3d')
	#ax.view_init(0, 0)
	import numpy
	def orthogonal_proj(zfront, zback):
	    a = (zfront+zback)/(zfront-zback)
	    b = -2*(zfront*zback)/(zfront-zback)
	    return numpy.array([[1,0,0,0],
	                        [0,1,0,0],
	                        [0,0,a,b],
	                        [0,0,0,zback]])
	proj3d.persp_transformation = orthogonal_proj
	ax.scatter(x,y,z)
	pl.savefig(filename,bbox_inches='tight',transparent=True) 
	pl.close()