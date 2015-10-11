import pandas as pd
df=pd.read_csv('a.csv',names=['date','value'])
df['date']=pd.to_datetime(df['date'])
dv=pd.DataFrame({'date':pd.date_range('20050104','20120228')})
df.reindex(index=['date'])
xx=pd.merge(dv,df,how='left',on='date') 
#xx['value']=xx['value'].astype('string')
xx['value']=xx['value'].fillna('.')
xx.to_csv('cc.csv',Index=None)
