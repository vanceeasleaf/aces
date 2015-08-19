import numpy as np
import pandas as pd
from aces.graph import plot,series
def getDos(frequencies,weights,filter=None):
	if filter is None:
		asum=np.array(np.cumsum(weights))
		filter=asum/asum[-1]<0.99
	freq=frequencies[filter]
	wei=weights[filter]
	frequency_points,sigma=get_draw_area(freq)
	#dos = np.array([get_density_of_states_at_freq(f,freq,wei,sigma)
	#					for f in frequency_points])
	dos=get_density_of_states_all(frequency_points,freq,wei,sigma)
	return frequency_points,dos
def getFilter(weights):
	asum=np.cumsum(weights,axis=1)
	filter=np.zeros_like(weights[0],dtype=bool)
	for a in asum:
		filter=a/a[-1]<0.99+filter

	return filter
def get_density_of_states_all(frequency_points,frequencies,weights,sigma):
	f0,f1=np.meshgrid(frequency_points,frequencies)
	ff=f0-f1
	peak=calc(sigma, ff)
	return np.dot(weights,peak)/np.sum(weights)
def get_density_of_states_at_freq( f,frequencies,weights,sigma):
	return np.sum(np.dot(
		weights, calc(sigma,frequencies - f))
		) /  np.sum(weights)
def calc(sigma, x):
	return 1.0 / np.sqrt(2 * np.pi) / sigma * \
		np.exp(-x**2 / 2.0 / sigma**2)		
def get_draw_area(frequencies,
				  freq_min=None,
				  freq_max=None,
				  freq_pitch=None):

	f_min = frequencies.min()
	f_max = frequencies.max()


	sigma = (f_max - f_min) /200
		
	if freq_min is None:
		f_min -= sigma * 10
	else:
		f_min = freq_min
		
	if freq_max is None:
		f_max +=sigma * 10
	else:
		f_max = freq_max

	if freq_pitch is None:
		f_delta = (f_max - f_min) / 1000
	else:
		f_delta = freq_pitch
	frequency_points = np.arange(f_min, f_max + f_delta * 0.1, f_delta)
	return frequency_points,sigma
def plot_smooth():		
	df=pd.read_csv("VDOS.txt",sep=r"[ \t]",engine="python");
	freq,ave=getDos(df['Freq_THz'],df['vdos_av'])
	dp=pd.DataFrame()
	dp['Freq_THz']=freq
	dp["vdos_av"]=ave
	dp.to_csv('smooth_dos.txt',sep='\t',index=False,float_format ="%f" );
	plot((freq,'Frequency(THz)'),(ave,'Phonon Density of states'),"smooth_dos.png",grid=True)	

	alldos=[df['vdos_x'],df['vdos_y'],df['vdos_z']]
	filter=getFilter(alldos)
	dp=pd.DataFrame()
	freq,ave=getDos(df['Freq_THz'],df['vdos_x'],filter)
	dp['Freq_x']=freq
	dp["vdos_x"]=ave
	freq,ave=getDos(df['Freq_THz'],df['vdos_y'],filter)
	dp['Freq_y']=freq
	dp["vdos_y"]=ave
	freq,ave=getDos(df['Freq_THz'],df['vdos_z'],filter)
	dp['Freq_z']=freq
	dp["vdos_z"]=ave
	dp.to_csv('smooth_dosxyz.txt',sep='\t',index=False,float_format ="%f" );	
	series(xlabel='Frequency(THz)',
		ylabel='Phonon Density of States',
		datas=[(dp['Freq_x'],dp["vdos_x"],"dos_x"),
		(dp['Freq_y'],dp["vdos_y"],"dos_y"),
		(dp['Freq_z'],dp["vdos_z"],"dos_z")]
		,filename='smooth_dosxyz.png',grid=True)
def plot_dos():
	df=pd.read_csv("VDOS.txt",sep=r"[ \t]",engine="python");
	freq,dosx,dosy,dosz=df['Freq_THz'],df['vdos_x'],df['vdos_y'],df['vdos_z']
	series(xlabel='Frequency(THz)',
		ylabel='Phonon Density of States',
		datas=[(freq,dosx,"dos_x"),
		(freq,dosy,"dos_y"),
		(freq,dosz,"dos_z")]
		,linewidth=0.3
		,filename='VDOS.png')
def select(x,n=1000):
	N=len(x)
	if N<n:
		return range(0,N)
	else:
		return range(0,N,N/n)
def plot_vacf():
	df=pd.read_csv("VACF.txt",sep=r"[ \t]",engine="python");
	time,vx,vy,vz=df['correlation_time(ps)'],df['vcaf_x'],df['vcaf_y'],df['vcaf_z']
	xx=select(time)
	time=time[xx]
	vx,vy,vz=vx[xx],vy[xx],vz[xx]
	series(xlabel='Correlation Time (ps)',
		ylabel='Normalized Velocity Auto Correlation Function',
		datas=[(time,vx,"vcf_x"),
		(time,vy,"vcf_y"),
		(time,vz,"vcf_z")]
		,linewidth=0.3
		,filename='VACF.png')

def plot_atomdos(atoms=range(24),filename='atom_dos.png'):
	import h5py
	db=h5py.File('velocity.hdf5')
	freq=db['/freq']
	datas=[]
	for i in atoms:
		dos=db['/dos_atom/%d'%i]
		x,y=getDos(freq,dos)
		datas.append((x,y,"dos_atom_%d"%i))
	series(xlabel='Frequency (THz)',
		ylabel='Phonon Density of States',
		datas=datas
		,linewidth=0.3
		,filename=filename,legend=False)
def plot_regiondos(regions,filename='region_dos.png'):
	import h5py
	db=h5py.File('velocity.hdf5')
	freq=db['/freq']
	datas=[]
	alldos=[]

	for region,rname in regions:
		dos=get_regiondos(region)
		alldos.append(dos)
	filter=getFilter(alldos)
	dp=pd.DataFrame()
	for region,rname in regions:
		dos=get_regiondos(region)
		x,y=getDos(freq,dos,filter)
		datas.append((x,y,"region:%s"%rname))
		dp['freq_%s'%rname]=x
		dp['dos_%s'%rname]=y
	dp.to_csv('region_dos.txt',sep='\t',index=False,float_format ="%f" );	
	series(xlabel='Frequency (THz)',
		ylabel='Phonon Density of States',
		datas=datas
		,linewidth=1
		,filename=filename,legend=True)
def get_regiondos(region):
	import h5py
	db=h5py.File('velocity.hdf5')
	dos=db['/dos_atom/%d'%region[0]]
	totaldos=np.zeros_like(dos)
	for i in region:
		totaldos+=db['/dos_atom/%d'%i]
	return totaldos
if __name__=='__main__':
	import json
	fp=open('qloops.txt')
	dirs=[]
	for line in fp:
		obj=json.loads(line)
		#if not obj['id']%3==0:continue
		dirs.append(dict(id=obj['id'],latysmall=obj['latysmall']))

	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from aces import tools
	ys=[]
	dp=pd.DataFrame()
	for dir in dirs:

		tools.cd(str(dir['id']))
		df=pd.read_csv("VDOS.txt",sep=r"[ \t]",engine="python");
		x,y=getDos(df['Freq_THz'],df['vdos_av'])
		dp['freq']=x;
		dp["latysmall=%s"%dir['latysmall']]=y
		tools.cd('..')
	dp.to_csv('total_dos.txt',sep='\t',index=False,float_format ="%f" );
	for dir in dirs:
		label="latysmall=%s"%dir['latysmall'];
		plt.plot(dp['freq'], dp[label],linewidth=1,label=label)
	plt.grid(True)
	plt.xlim(min(x), max(x))
	plt.xlabel('Frequency(THz)')
	plt.ylabel('Density of states')
	plt.legend(loc='best').get_frame().set_alpha(0.0)
	plt.savefig('smoothVDOS.png',bbox_inches='tight',transparent=True)