import pandas as pd
df=pd.read_csv("region_dos.txt",sep=r"[ \t]",engine="python");
npair=len(df.columns)/2
datas=[]
for i in range(npair):
	rname=df.columns[i*2][5:]
	datas.append((df['freq_'+rname],df['dos_'+rname],"region:"+rname))
dc=pd.read_csv("graphenedos.txt",sep=r"[ \t]",engine="python");
datas.append((dc[dc.columns[0]],dc[dc.columns[1]],'GNR'))
from aces.graph import series
series(xlabel='Frequency (THz)',
	ylabel='Phonon Density of States',
	datas=datas
	,linewidth=1
	,filename='camparedos.png',legend=True,xmax=60)