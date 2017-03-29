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
if False:
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
		for var in ['tasmax','tasmin']:
			try:
				nc_hist=Dataset(files_to_treat[model][var+'_historical'])
				datevar = []
				datevar.append(num2date(nc_hist.variables['time'][:],units = nc_hist.variables['time'].units,calendar = nc_hist.variables['time'].calendar))
				year_hist=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])

				nc_45=Dataset(files_to_treat[model][var+'_rcp45'])
				datevar = []
				datevar.append(num2date(nc_45.variables['time'][:],units = nc_45.variables['time'].units,calendar = nc_45.variables['time'].calendar))
				year_45=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])

				files_to_treat[model][var+'_year']=np.concatenate((year_hist,year_45),axis=0)

				files_to_treat[model]['lat']=nc_hist.variables['lat'][:]
				files_to_treat[model]['lon']=nc_hist.variables['lon'][:]

				var_hist=nc_hist.variables[var][:,:,:]
				var_45=nc_45.variables[var][:,:,:]
				var_in=np.concatenate((var_hist,var_45),axis=0)
				if var_in.mean()>150:var_in-=273.15
				files_to_treat[model][var+'_values']=var_in
			except:
				print model



	for model in files_to_treat.keys():
		GMT=glob.glob('../wlcalculator/data/cmip5_ver002/'+model.lower()+'.rcp45.r1i*')
		if len(GMT)>0:
			files_to_treat[model]['GMT']=GMT[0]

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
				#print change,np.nanmin(abs(rmean-change)),year[closest],rmean[closest]
				if np.nanmin(abs(rmean-change))<0.1:
					files_to_treat[model][change]=year[closest]

			files_to_treat[model]['GMT_year']=year
			files_to_treat[model]['GMT_values']=GMT
		except:
			print model

	for model in files_to_treat.keys():
		if len(files_to_treat[model].keys())<16:
			files_to_treat.pop(model, None)


	# median
	for var in ['tasmax','tasmin']:
		min_yr=[]
		max_yr=[]
		for model in files_to_treat.keys():
			min_yr.append(files_to_treat[model][var+'_year'].min())
			max_yr.append(files_to_treat[model][var+'_year'].max())
		print var, min_yr, max_yr

	for qu in [25,50,75]:
		files_to_treat[qu]={}

		for var in ['tasmax','tasmin']:
			for model in files_to_treat.keys():
				if model not in [25,50,75]:			
					os.system('cdo -O mergetime '+files_to_treat[model][var+'_historical']+' '+files_to_treat[model][var+'_rcp45']+' CMIP5_regrid_ensemble_selection/'+var+'_'+model+'_73x96.nc')
					files_to_treat[model][var+'_merged']='CMIP5_regrid_ensemble_selection/'+var+'_'+model+'_73x96.nc'

			command='cdo -O enspctl,'+str(qu)+' '
			for model in files_to_treat.keys():
				if model not in [25,50,75]:	
					os.system('cdo -O delete,timestep=-1 '+files_to_treat[model][var+'_merged']+' CMIP5_regrid_ensemble_selection/'+var+'_'+model+'_73x96_cut.nc')
					if model in ['HadGEM2-CC']:
						os.system('rm CMIP5_regrid_ensemble_selection/'+var+'_'+model+'_73x96_cut.nc ')
						os.system('cp CMIP5_regrid_ensemble_selection/'+var+'_'+model+'_73x96.nc CMIP5_regrid_ensemble_selection/'+var+'_'+model+'_73x96_cut.nc ')
					command+='CMIP5_regrid_ensemble_selection/'+var+'_'+model+'_73x96_cut.nc '
			command+='CMIP5_regrid_ensemble_selection/'+var+'_'+str(qu)+'_73x96_cut.nc'
			os.system(command)
			files_to_treat[qu][var+'_merged']='CMIP5_regrid_ensemble_selection/'+var+'_'+str(qu)+'_73x96_cut.nc'

	# median GMT
	min_yr=[]
	max_yr=[]
	for model in files_to_treat.keys():
		if model not in [25,50,75]:
			min_yr.append(files_to_treat[model]['GMT_year'].min())
			max_yr.append(files_to_treat[model]['GMT_year'].max())	
	print min_yr, max_yr

	time_axis=range(int(max(min_yr)),int(min(max_yr)),1)


	# was mean! not median change
	GMTs=np.zeros([len(time_axis),13])*np.nan	
	count=0
	for model in files_to_treat.keys():
		if model not in [25,50,75]:
			relevant_years=np.where((files_to_treat[model]['GMT_year']>=time_axis[0]) & (files_to_treat[model]['GMT_year']<=time_axis[-1]))[0]
			if model in ['HadGEM2-ES','HadGEM2-CC']:
				relevant_before_2005=np.where((files_to_treat[model]['GMT_year']>=time_axis[0]) & (files_to_treat[model]['GMT_year']<=2005))[0]
				relevant_after_2006=np.where((files_to_treat[model]['GMT_year']>=2006) & (files_to_treat[model]['GMT_year']<=time_axis[-1]))[0]
				yr_2006=np.mean([files_to_treat[model]['GMT_values'][relevant_before_2005][-1],files_to_treat[model]['GMT_values'][relevant_after_2006][0]])
				GMT_fixed=np.array(list(files_to_treat[model]['GMT_values'][relevant_before_2005])+[yr_2006]+list(files_to_treat[model]['GMT_values'][relevant_after_2006]))
				GMTs[:,count]=GMT_fixed
				print GMT_fixed
			else:	
				print model
				GMTs[:,count]=files_to_treat[model]['GMT_values'][relevant_years]
			count+=1


	for qu in [25,50,75]:
		# time slice
		ave_window=20
		GMT=np.percentile(GMTs,qu,axis=1)
		rmean=GMT.copy()*np.nan
		for i in range(19,len(rmean)):
			#print i-ave_window+1,i
			rmean[i]=GMT[i-ave_window:i].mean()
		#rmean_dict[ds]=rmean-rmean.ix[2015].values

		rmean=rmean-rmean[np.where(year==2010)[0]]

		for change in [-0.5,+0.811,+1.311]:
			closest=np.nanargmin(abs(rmean-change))
			#print change,np.nanmin(abs(rmean-change)),year[closest],rmean[closest]
			if np.nanmin(abs(rmean-change))<0.1:
				files_to_treat[qu][change]=year[closest]





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


files_to_treat.pop(25, None)
files_to_treat.pop(50, None)
files_to_treat.pop(75, None)

cmip5_dict={}
TXx_0p5_hist=np.zeros([94,14])*np.nan

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
        cmip5_dict[model][var].derive_regional_masking(region_type='region',dataset='HadEX2',overwrite=False)

        #cmip5_dict[model][var]._masks['global']=cmip5_dict['HadEX2'][var]._masks['global']
        
        # Derive time slices
        cmip5_dict[model][var].derive_time_slices(ref_period,target_periods,period_names)
        cmip5_dict[model][var].derive_distributions()

        cmip5_dict[model][var].derive_pdf_difference(str(-0.5),'ref',pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])

        if model in [25,50,75]:linewidth=3
    	else: linewidth=1
        plt.plot(cmip5_dict[model][var]._distributions['global']['pdf']['xaxis'],cmip5_dict[model][var]._distributions['global']['pdf']['ref_-0.5'],linewidth=linewidth)
        


        TXx_0p5_hist[:,mod_index]=cmip5_dict[model][var]._distributions['global']['diff']['ref_-0.5']

        #cmip5_dict[model][var].derive_pdf_difference(str(+0.811),str(+1.311),pdf_method=pdf_method,bin_range=varin_dict[var]['cut_interval'])
        #plt.plot(cmip5_dict[model][var]._distributions['Europe']['pdf']['xaxis'],cmip5_dict[model][var]._distributions['Europe']['pdf']['1.311_0.811'])

       	#plt.show()
        #asdasd

plt.show()



