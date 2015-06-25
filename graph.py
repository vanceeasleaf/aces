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
	pl.close()
def plot(x,y,filename,grid=False):
	series(x[1],y[1],[(x[0],y[0],y[1])],filename,legend=False,grid=grid)
	
def series(xlabel,ylabel,datas,filename,linewidth=1,legend=True,grid=False):
	pl.figure()
	pl.xlabel(xlabel)
	pl.ylabel(ylabel)
	if len(datas)==1:
		serie=datas[0]
		pl.plot(serie[0],serie[1],label=serie[2],linewidth=linewidth,color='r')
	else:
		for serie in datas:
			pl.plot(serie[0],serie[1],label=serie[2],linewidth=linewidth)
	if legend:
		pl.legend(loc='best').get_frame().set_alpha(0.0)
	if grid:
		pl.grid(True)
	pl.savefig(filename,bbox_inches='tight',transparent=True) 
	pl.close()

def surf(X, Y, Z,filename):
	from mpl_toolkits.mplot3d import Axes3D
	fig = pl.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X, Y, Z)
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