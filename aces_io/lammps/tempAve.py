
from  aces.graph import twinx
import pandas as pd
def drawTempAve():
	df=pd.read_csv('tempAve.txt',sep='\t')
	twinx((df['Coord'],'x(Augstrom)'),(df['Temp'],'temperature(K)'),(df['Jx'],'heat flux'),'tempAve.png')
 
if __name__=='__main__':
	drawTempAve()