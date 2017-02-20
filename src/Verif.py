import numpy as np
import pickle
import datetime
from netCDF4 import Dataset
import pandas as pd
import scipy.io

SAVEPATH='tmp/'

path='/mnt/web/brey/hincast/'

SAVEPATH='/mnt/pandora/Users_Critech/Brey/'

obsname='/mnt/pandora/Users_Critech/Brey/Croton.mat'

basename='med'

t1='20131201.00'

t2='20160229.12'

sdate=datetime.datetime.strptime(t1,'%Y%m%d.%H')

edate=datetime.datetime.strptime(t2,'%Y%m%d.%H')

path0=path+'{}/'.format(t1)

with open(path0+basename+'.pkl', 'r') as f:
    ptr=pickle.load(f)

iloc=1844 # CHANGE THIS FOR ANALYSIS OF ANOTHER OBS POINT
location='Crotone'

# ## Get All forecasting data

dt=(edate-sdate).total_seconds()
ndt=dt/(3600*12)
ndt=np.int(ndt)+1

combined=[] # store time and height
tstamp = [] # for use later

tot=[]

indx=pd.date_range(sdate, edate, freq='12H')


val=[]

for it in range(ndt):
    idate=sdate+datetime.timedelta(hours=12*it)
    dstamp=datetime.datetime.strftime(idate,'%Y%m%d.%H')
    path0=path+'{}/'.format(dstamp)
    filename=path0+'trih-'+basename+'.nc'
    d =  Dataset(filename)
    ha=d.variables['ZWL'][:,ptr[iloc]] # all values
    t=d.variables['time'][:]
    tw=[]
    for it in t:
        tw.append(idate+datetime.timedelta(seconds=np.int(it)))
    ttw=[(item-t[0])/60. for item in t]
    dic={'time':tw, 'numerical storm surge':ha, 'forecast_time':ttw}
    val.append(ha)
    data=pd.DataFrame.from_dict(dic)
    data=data.set_index('time')
    tot.append(data)

tota=pd.concat(tot,keys=indx)

tota.drop('forecast_time',1, inplace=True)

tota.index.tolist()[:10]


ht=np.array([np.arange(4321) for i in range(tota.index.levels[0].shape[0])]).ravel()

len(tota.index.tolist())

index0=[(a, c) for (a,b),c in zip(tota.index.tolist(),ht)]

tota.index=pd.MultiIndex.from_tuples(index0,names=tota.index.names)

s=tota.unstack(level=0)

# #### output for verif

# ### rolling average for all

av5=s.rolling(5,center=True, min_periods=3).mean()[::60][:25] # mean within 5 min for every hour centered

av5.index=av5.index/60

sav5=av5.unstack()

sav5.index=sav5.index.droplevel(level=0)

#s24=av5.iloc[:,::2].unstack()

#s24.index=s24.index.droplevel(level=0)

# #### read cleaned obs

mat=scipy.io.loadmat(obsname)


obs=pd.DataFrame({'time':mat['times'].flatten(),'s':mat['ss_Croton'].flatten(), 'offset':mat['offset'].flatten()})

obs['time']=pd.to_datetime(obs['time'].str[:-4])

#take out the duplicate 24h value
obs=obs[obs.offset != 24]

#dublicate for x.12 analysis
obs12s=obs['s'][12:]
obs12times=obs['time']+datetime.timedelta(hours=12)
obs12offset=obs['offset']
obs12=pd.DataFrame({'time':obs12times.values[:-12],'s':obs12s.values,'offset':obs12offset.values[:-12]})

#add the 24 forecast for the 00

#copy the 0 values
ext=obs[::24].copy()

#change the offset
ext.offset=24

#change time to the previous day
ext.time=ext.time-datetime.timedelta(days=1)

ext=ext.drop(0) # drop the first value which is not valid

ext=ext.set_index(['time','offset'])

obs=obs.set_index(['time','offset'])

obs_ex=pd.concat([obs,ext]).sort_index(level=0)

#add the 24 for the 12
#copy the 0 values
ext12=obs12[::24].copy()

#change the offset
ext12.offset=24

#change time to the previous day
ext12.time=ext12.time-datetime.timedelta(days=1)

ext12=ext12.drop(0)

ext12=ext12.set_index(['time','offset'])

obs12=obs12.set_index(['time','offset'])

obs12_ex=pd.concat([obs12,ext12]).sort_index(level=0)


#join them
obsf=pd.concat([obs_ex,obs12_ex]).sort_index(level=0)

obsf['s']=obsf['s'].fillna('NaN')


# ### join the 2 dataframes

ver=pd.concat([obsf,sav5],axis=1)

output=ver.loc['2013-12-01 00:00:00':'2016-02-29 12:00:00']


# #### output

with open(SAVEPATH+'/{}.csv'.format(location), 'w') as f:
    f.write('# variable: s\n')
    f.write('# units: $m$\n')
    f.write('date offset obs fcst\n')

with open(SAVEPATH+'{}.csv'.format(location), 'a') as f:
        output.to_csv(f, date_format='%Y%m%d%H', sep='\t', header=None)


