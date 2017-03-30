


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
#plt.rcParams['figure.figsize'] = 4,6
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


with open('../CMIP5/varoutdict_cmip5_hadex2-grid.pkl', 'rb') as input:
   cmip5_dict = pickle.load(input)   



time_slice_dict={
  1:{'name':'ref_-0.5','longname':'ref vs -0.5'},
  2:{'name':'0.811_0.311','longname':'+1.5 vs +1'},
  3:{'name':'1.311_0.811','longname':'+2 vs +1.5'},
  4:{'name':'1.811_1.311','longname':'+2.5 vs +2'}
}

f,pl=plt.subplots(nrows=4,ncols=2,figsize=(10,12))
pplot=pl.flatten()
count=0
for ID in time_slice_dict.keys():
	for var in ['TXx','TNn']:
		PDFs=np.zeros([512,N_model])*np.nan
		for model,mod_index in zip(files_to_treat.keys(),range(N_model)):
			try:
				pplot[count].plot(cmip5_dict[model][var]._distributions['global']['pdf']['xaxis'],cmip5_dict[model][var]._distributions['global']['pdf'][time_slice_dict[ID]['name']],linewidth=0.3,color='gray')
				PDFs[:,mod_index]=cmip5_dict[model][var]._distributions['global']['pdf'][time_slice]
			except:
				print time_slice

		pl95=np.nanpercentile(PDFs,95,axis=1)
		pl5=np.nanpercentile(PDFs,5,axis=1)
		pplot[count].fill_between(cmip5_dict[model][var]._distributions['global']['pdf']['xaxis'],
		              pl95,pl5,color='gray',
		                    alpha=0.25)

		pplot[count].plot([0,0],[0,1],color='k',linestyle=':')		
		pplot[count].set_ylim((0,0.013))

		if (-1)**count>0:
			pplot[count].set_ylabel(time_slice_dict[ID]['longname'])
			pplot[count].set_yticks([], [])
		if (-1)**count<0:
			pplot[count].yaxis.tick_right()
		if count<2: pplot[count].set_title(varin_dict[var]['longname'])
		if count>5: 
			pplot[count].set_xlabel(varin_dict[var]['unit'])
		if count<6:
			pplot[count].set_xticks([], [])
		count+=1

plt.subplots_adjust(wspace=0,hspace=0)
plt.savefig('../CMIP5/CMIP5_time_slices.png',dpi=300)
plt.clf()






