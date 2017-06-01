datas=[]
import pandas as pd
import numpy as np
from aces.tools import *
youngs=[]
temps=range(200,300,10)
for i in range(10):
	print i
	cd(str(i))
	from aces.App import App
	App().runner.post()
	
	df=pd.read_csv("cal_stress.txt",sep=r"[ \t]",engine="python");
	datas.append([df['Strain'],df['Stress_GPa'],str(temps[i])+'K'	])
	youngs.append(float(read("YoungsModulus.txt")))
	cd('..')
from aces.graph import series,plot
plot((np.array(temps),'Temperature (K)'),(np.array(youngs),'Young\'s Modulus (GPa)'),'youngs.png')
series(xlabel='Strain',
	ylabel='Stress (GPa)',
	datas=datas
	,linewidth=1
	,filename='stress_strain.png',legend=True)