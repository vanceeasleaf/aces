import numpy as np
import pandas as pd
from aces.graph import plot
def getDos(frequencies,weights):
	filter=weights>weights.max()/10000
	freq_max=frequencies[filter].max()
	frequency_points,sigma=get_draw_area(frequencies,freq_max=freq_max)
	dos = np.array([get_density_of_states_at_freq(f,frequencies,weights,sigma)
						for f in frequency_points])
	return frequency_points,dos
	
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


	sigma = (f_max - f_min) /100
		
	if freq_min is None:
		f_min -= sigma * 10
	else:
		f_min = freq_min
		
	if freq_max is None:
		f_max +=sigma * 10
	else:
		f_max = freq_max

	if freq_pitch is None:
		f_delta = (f_max - f_min) / 200
	else:
		f_delta = freq_pitch
	frequency_points = np.arange(f_min, f_max + f_delta * 0.1, f_delta)
	return frequency_points,sigma
def plot_smooth():		
	df=pd.read_csv("VDOS.txt",sep=r"[ \t]",engine="python");
	x,y=getDos(df['Freq_THz'],df['vdos_av'])
	dp=pd.DataFrame()
	dp['freq']=x;
	dp["dos"]=y
	dp.to_csv('smooth_dos.dat',sep='\t',index=False,float_format ="%f" );
	plot((x,'Frequency(THz)'),(y,'Density of states'),"smooth_dos.png",grid=True)	
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