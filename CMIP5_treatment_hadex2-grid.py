import numpy as np
import matplotlib.pyplot as plt
import dimarray as da 
from netCDF4 import Dataset,netcdftime,num2date
import sys,glob,datetime,pickle,os,itertools



os.chdir('/Users/peterpfleiderer/Documents/Projects/0p5_observed/')

# #regrid
# grid='73x96'

# mygrid=open('CMIP5/CMIP5_regrid_HadEX2/'+grid+'.txt','w')
# mygrid.write('gridtype=lonlat\nxsize='+str(96)+'\nysize='+str(73)+'\nxfirst='+str(0.0)+'\nxinc='+str(3.75)+'\nyfirst='+str(-90)+'\nyinc='+str(2.5))
# mygrid.close()

# all_files=glob.glob('cmip5_Xev_from_Erich_Fischer/*')
# for file in all_files:
# 	os.system('cdo remapbil,CMIP5/CMIP5_regrid_HadEX2/'+grid+'.txt '+file+' CMIP5/CMIP5_regrid_HadEX2/'+file.split('/')[-1].replace('.nc','_'+grid+'.nc'))	


# only few models, could be improved!
if True:
	files_to_treat={}
	model_dict={}
	for var in ['tasmax','tasmin']:
		for scenario in ['rcp45','historical']:
			all_files = glob.glob('CMIP5/CMIP5_regrid_HadEX2/'+var+'*'+scenario+'*r1i*')
			for file in all_files:
				model=file.split('/')[-1].split('_')[1]
				if model not in files_to_treat.keys():
					files_to_treat[model]={}
				if var not in files_to_treat[model].keys():
					files_to_treat[model][var+'_'+scenario]=file

		for model in files_to_treat.keys():
			try:
				os.system('cdo -O mergetime '+files_to_treat[model][var+'_historical']+' '+files_to_treat[model][var+'_rcp45']+' CMIP5/CMIP5_regrid_HadEX2_ensemble_selection/'+var+'_'+model+'_73x96.nc')
				files_to_treat[model][var+'_merged']='CMIP5/CMIP5_regrid_HadEX2_ensemble_selection/'+var+'_'+model+'_73x96.nc'
			except: 
				pass

	for model in files_to_treat.keys():
		model_dict[model]=[]
		GMT=glob.glob('../wlcalculator/data/cmip5_ver002/'+model.lower()+'.rcp45.r1i*')
		if len(GMT)>0:
			files_to_treat[model]['GMT']=GMT[0]

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
			try:
				ref_warming=rmean[np.where(year==2010)[0]][0]
			except:
				ref_warming=rmean[np.where(year==2010)[0]]
			rmean=rmean-ref_warming

			for change in [-0.5,+0.311,+0.811,+1.311,+1.811]: 
				closest=np.nanargmin(abs(rmean-change))
				print change,np.nanmin(abs(rmean-change)),year[closest],rmean[closest]
				if np.nanmin(abs(rmean-change))<0.1:
					files_to_treat[model][change]=year[closest]
				else:
					model_dict[model]+=[change]

		else:
			model_dict[model]+=['no GMT']



	for model in files_to_treat.keys(): 
		if len(files_to_treat[model].keys())<10:
			model_dict[model]+=['missing input']
			files_to_treat.pop(model, None)


	with open('CMIP5/cmip5_HadEX2_time_slices_and_files.pkl', 'wb') as output:
	    pickle.dump(files_to_treat, output, pickle.HIGHEST_PROTOCOL)

data_missing=''
ensemble=''
for model in model_dict.keys():
	if 'missing input' in model_dict[model]:
		data_missing+=model+' '
	else:
		ensemble+=model+' '


adsasd

# IMPORT AND CONFIG 
import numpy as np
import netCDF4 as net
import dimarray as da 
import sys,glob,datetime,pickle,os,itertools
import pandas as pd
import matplotlib.pylab as plt 
from netCDF4 import Dataset,netcdftime,num2date
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon


#plt.style.use('ggplot')
#plt.rcParams['figure.figsize'] = 8,6
from matplotlib import rc
rc('text', usetex=True)

os.chdir('/Users/peterpfleiderer/Documents/Projects/0p5_observed/pdf_processing/')
try:
    import pdf_processing as pdf; reload(pdf)
except ImportError:
    raise ImportError(
        "cannot find PDF_Processing code")

###########
# Settings
###########

# PDF Method (currently defined: hist, python_silverman)
pdf_method='python_silverman'

# variables
varin_dict={
    'TXx':{'var_name':'tasmax','longname':'Hot extremes (TXx)','unit':'TXx [$^\circ$ C]','cut_interval':[-2,3]},
    'TNn':{'var_name':'tasmin','longname':'Cold extremes (TNn)','unit':'TNn [$^\circ$ C]','cut_interval':[-3,5]}
    }

N_model=len(files_to_treat.keys())


if True:
	cmip5_dict={}

	for var in ['TXx','TNn']:
		for model,mod_index in zip(files_to_treat.keys(),range(N_model)):
			print model
			if model not in cmip5_dict.keys(): cmip5_dict[model]={}

			# time informations and periods
			ref_period=[1991,2010]
			target_periods=[]
			period_names=[]
			for change in [-0.5,+0.311,+0.811,+1.311,+1.811]:
				try:
					end_year=files_to_treat[model][change]
					target_periods.append([end_year-19,end_year])
					period_names.append(str(change))
				except:
					print change, 'not reached'
			target_periods.append(ref_period)
			period_names.append('ref')

			os.chdir('/Users/peterpfleiderer/Documents/Projects/0p5_observed/')
			nc_in=Dataset(files_to_treat[model][varin_dict[var]['var_name']+'_merged'])
			datevar = []
			datevar.append(num2date(nc_in.variables['time'][:],units = nc_in.variables['time'].units,calendar = nc_in.variables['time'].calendar))
			year=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])

			lat=nc_in.variables['lat'][:]
			lon=nc_in.variables['lon'][:]

			# combine datasets
			var_name=varin_dict[var]['var_name']
			var_in=nc_in.variables[var_name][:,:,:]
			if var_in.mean()>150:var_in-=273.15

			input_data=da.DimArray(var_in[:,:,:].squeeze(), axes=[year, lat, lon],dims=['year', 'lat', 'lon'] )

			# Mask for data availability (individual for each model)
			os.chdir('/Users/peterpfleiderer/Documents/Projects/0p5_observed/pdf_processing/')
			cmip5_dict[model][var]=pdf.PDF_Processing(var)
			cmip5_dict[model][var].mask_for_ref_period_data_coverage(input_data,ref_period,check_ref_period_only=False,target_periods=target_periods,dataset='HadEX2')
			cmip5_dict[model][var].derive_regional_masking(region_type='EU_NA_AS_RU_AU',dataset='HadEX2',overwrite=False)

			# Derive time slices
			cmip5_dict[model][var].derive_time_slices(ref_period,target_periods,period_names)
			cmip5_dict[model][var].derive_distributions()

			cmip5_dict[model][var].derive_pdf_difference(str(-0.5),'ref',pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
			if str(0.311) in period_names:
				cmip5_dict[model][var].derive_pdf_difference('ref',str(0.311),pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
				if str(0.811) in period_names:
					cmip5_dict[model][var].derive_pdf_difference(str(0.311),str(0.811),pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
					if str(1.311) in period_names:
						cmip5_dict[model][var].derive_pdf_difference(str(0.811),str(1.311),pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
						if str(1.811) in period_names:
							try:
								cmip5_dict[model][var].derive_pdf_difference(str(1.311),str(1.811),pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
							except:
								pass
				

		
	with open('../CMIP5/varoutdict_cmip5_hadex2-grid.pkl', 'wb') as output:
	    pickle.dump(cmip5_dict, output, pickle.HIGHEST_PROTOCOL)
	    
    





