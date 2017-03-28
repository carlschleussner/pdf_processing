import numpy as np
import matplotlib.pyplot as plt
import dimarray as da 
from netCDF4 import Dataset,netcdftime,num2date
import sys,glob,datetime,pickle,os,itertools



os.chdir('/Users/peterpfleiderer/Documents/Projects/0p5_observed/')

# #regrid
# grid='73x96'

# mygrid=open('CMIP5_regrid/'+grid+'.txt','w')
# mygrid.write('gridtype=lonlat\nxsize='+str(96)+'\nysize='+str(73)+'\nxfirst='+str(0.0)+'\nxinc='+str(3.75)+'\nyfirst='+str(-90)+'\nyinc='+str(2.5))
# mygrid.close()

# all_files=glob.glob('cmip5_Xev_from_Erich_Fischer/*')
# for file in all_files:
# 	os.system('cdo remapbil,CMIP5_regrid/'+grid+'.txt '+file+' CMIP5_regrid/'+file.split('/')[-1].replace('.nc','_'+grid+'.nc'))	




files_to_treat={}

for var in ['tasmax','tasmin']:
	for scenario in ['rcp45','historical']:
		all_files = glob.glob('CMIP5_regrid/'+var+'*'+scenario+'*r1i*')
		for file in all_files:
			model=file.split('/')[-1].split('_')[1]
			if model not in files_to_treat.keys():
				files_to_treat[model]={}
			if var not in files_to_treat[model].keys():
				files_to_treat[model][var+'_'+scenario]=file


for model in files_to_treat.keys():
	GMT=glob.glob('../wlcalculator/data/cmip5_ver002/'+model.lower()+'.rcp45.r1i*')
	if len(GMT)>0:
		files_to_treat[model]['GMT']=GMT[0]
	if len(files_to_treat[model].keys())<5:
		files_to_treat.pop(model, None)


for model in files_to_treat.keys():
	try:
		nc_in=Dataset(files_to_treat[model]['GMT'],"r")

		# handle time information
		time=nc_in.variables['time'][:]
		year=time/12

		# GMT
		GMT = nc_in.variables['tas_global'][:]	

		ave_window=20

		rmean=GMT.copy()*np.nan
		for i in range(19,len(rmean)):
			#print i-ave_window+1,i
			rmean[i]=GMT[i-ave_window:i].mean()
		#rmean_dict[ds]=rmean-rmean.ix[2015].values

		rmean=rmean-rmean[np.where(year==2010)[0]]

		for change in [-0.5,+0.811,+1.311]:
			closest=np.nanargmin(abs(rmean-change))
			if np.nanmin(abs(rmean+0.5))<0.1:
				files_to_treat[model][change]=year[closest]
	except:
		print model


