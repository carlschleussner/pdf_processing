


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

colordict={
    'ERA':"#49bc2f",
    '20CR':"#fb9b32"}

with open('../varoutdict_10000_reanalysis.pkl', 'rb') as input:
    varoutdict_rea = pickle.load(input)   

with open('../CMIP5/varoutdict_cmip5_era-grid.pkl', 'rb') as input:
   cmip5_dict = pickle.load(input)   





f=plt.figure(figsize=(4,4))
for var in ['TXx']:
	PDFs=np.zeros([512,N_model])*np.nan
	for model,mod_index in zip(files_to_treat.keys(),range(N_model)):
		plt.plot(cmip5_dict[model][var]._distributions['global']['pdf']['xaxis'],cmip5_dict[model][var]._distributions['global']['pdf']['ref_-0.5'],linewidth=0.3,color='gray')
		PDFs[:,mod_index]=cmip5_dict[model][var]._distributions['global']['pdf']['ref_-0.5']

	pl95=np.percentile(PDFs,95,axis=1)
	pl5=np.percentile(PDFs,5,axis=1)
	plt.fill_between(cmip5_dict[model][var]._distributions['global']['pdf']['xaxis'],
                        pl95,pl5,color='gray',
                              alpha=0.25)

	for dataset in ['ERA']:
		plt.plot(varoutdict_rea[dataset][var]._distributions['global']['pdf']['xaxis'],
                  varoutdict_rea[dataset][var]._distributions['global']['pdf']['Recent_ref'],
                  label=dataset,
                  color=colordict[dataset],
                  linewidth=2)


	plt.xlabel(varin_dict[var]['unit'])
	plt.savefig('../CMIP5/CMIP5_'+var+'_ERA.png',dpi=300)
	plt.clf()
#plt.show()


# # check if masks are correct
# fig,pl=plt.subplots(nrows=1,ncols=2,figsize=(8,3))
# pplot=pl.flatten()
# mask=np.array(varoutdict_rea['ERA']['TXx']._masks['global'].copy())
# mask[np.isfinite(mask)]=1
# varoutdict_rea['ERA']['TXx'].plot_map(mask,
#                                     ax=pplot[0],
#                                     color_bar=False,
#                                     show=False)
# pplot[1].set_title('ERA')

# example_model=cmip5_dict.keys()[0]
# mask=np.array(cmip5_dict[example_model]['TXx']._masks['global'].copy())
# mask[np.isfinite(mask)]=1
# cmip5_dict[example_model]['TXx'].plot_map(mask,
#                                     ax=pplot[1],
#                                     color_bar=False,
#                                     show=False)
# pplot[1].set_title(example_model)

# plt.tight_layout()
# plt.savefig('../CMIP5/mask_ERA_CMIP5.png',dpi=300)
# plt.clf()

