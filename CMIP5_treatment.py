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


# only 13 models, could be improved!
if True:
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

			print model,rmean

			for change in [-0.5,+0.811,+1.311]:
				closest=np.nanargmin(abs(rmean-change))
				print change,np.nanmin(abs(rmean-change)),year[closest],rmean[closest]
				if np.nanmin(abs(rmean-change))<0.1:
					files_to_treat[model][change]=year[closest]
		except:
			print model

	for model in files_to_treat.keys():
		if len(files_to_treat[model].keys())<8:
			files_to_treat.pop(model, None)

	# # model median
	# files_to_treat['median']={}
	# for scenario in ['rcp45','historical']:
	# 	TXx_command='cdo enspctl,50 '
	# 	TNn_command='cdo enspctl,50 '
	# 	GMT_command='cdo enspctl,50 '
	# 	for model in files_to_treat.keys():
	# 		if model not in ['median','HadGEM2-CC','CCSM4']:
	# 			TXx_command+=files_to_treat[model]['tasmax'+'_'+scenario]+' '
	# 			TNn_command+=files_to_treat[model]['tasmin'+'_'+scenario]+' '
	# 			GMT_command+=files_to_treat[model]['GMT']+' '
	# 			os.system('cdo info -seldate '+files_to_treat[model]['GMT'])

	# 	TXx_command+='CMIP5_regrid/tasmax_median_'+scenario+'_73x96.nc'
	# 	TNn_command+='CMIP5_regrid/tasmin_median_'+scenario+'_73x96.nc'
	# 	GMT_command+='CMIP5_regrid/GMT_median.nc'
	# 	os.system(TXx_command)
	# 	os.system(TNn_command)
	# 	os.system(GMT_command)



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


plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 8,6
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




cmip5_dict={}
TXx_0p5_hist=np.zeros([512,14])*np.nan

for model,mod_index in zip(files_to_treat.keys(),range(13)):
	print model
	cmip5_dict[model]={}

	# time informations and periods
	ref_period=[1991,2010]
	target_periods=[]
	period_names=[]
	for change in [-0.5,+0.811,+1.311]:
		end_year=files_to_treat[model][change]
		print end_year
		target_periods.append([end_year-19,end_year])
		period_names.append(str(change))
	target_periods.append(ref_period)
	period_names.append('ref')

	for var in ['TXx']:
		os.chdir('/Users/peterpfleiderer/Documents/Projects/0p5_observed/')
		nc_hist=Dataset(files_to_treat[model][varin_dict[var]['var_name']+'_historical'])
		datevar = []
		datevar.append(num2date(nc_hist.variables['time'][:],units = nc_hist.variables['time'].units,calendar = nc_hist.variables['time'].calendar))
		year_hist=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])

		nc_45=Dataset(files_to_treat[model][varin_dict[var]['var_name']+'_rcp45'])
		datevar = []
		datevar.append(num2date(nc_45.variables['time'][:],units = nc_45.variables['time'].units,calendar = nc_45.variables['time'].calendar))
		year_45=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])

		year=np.concatenate((year_hist,year_45),axis=0)

		lat=nc_hist.variables['lat'][:]
		lon=nc_hist.variables['lon'][:]

		# combine datasets
		var_name=varin_dict[var]['var_name']
		var_hist=nc_hist.variables[var_name][:,:,:]
		var_45=nc_45.variables[var_name][:,:,:]
		var_in=np.concatenate((var_hist,var_45),axis=0)
		if var_in.mean()>150:var_in-=273.15

		input_data=da.DimArray(var_in[:,:,:].squeeze(), axes=[year, lat, lon],dims=['year', 'lat', 'lon'] )

		# Mask for data availability (individual for each model)
		os.chdir('/Users/peterpfleiderer/Documents/Projects/0p5_observed/pdf_processing/')
		cmip5_dict[model][var]=pdf.PDF_Processing(var)
		cmip5_dict[model][var].mask_for_ref_period_data_coverage(input_data,ref_period,check_ref_period_only=False,target_periods=target_periods,dataset='HadEX2')
        cmip5_dict[model][var].derive_regional_masking(region_type='region',dataset='HadEX2',overwrite=False)

        #cmip5_dict[model][var]._masks['global']=cmip5_dict['HadEX2'][var]._masks['global']
        
        # Derive time slices
        cmip5_dict[model][var].derive_time_slices(ref_period,target_periods,period_names)
        cmip5_dict[model][var].derive_distributions()

        cmip5_dict[model][var].derive_pdf_difference(str(-0.5),'ref',pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
        #plt.plot(cmip5_dict[model][var]._distributions['Europe']['pdf']['xaxis'],cmip5_dict[model][var]._distributions['Europe']['pdf']['ref_-0.5'])
        TXx_0p5_hist[:,mod_index]=cmip5_dict[model][var]._distributions['Europe']['pdf']['ref_-0.5']

        #cmip5_dict[model][var].derive_pdf_difference(str(+0.811),str(+1.311),pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
        #plt.plot(cmip5_dict[model][var]._distributions['Europe']['pdf']['xaxis'],cmip5_dict[model][var]._distributions['Europe']['pdf']['1.311_0.811'])

       	#plt.show()
        #asdasd





