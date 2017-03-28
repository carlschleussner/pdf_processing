
import numpy as np
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import datetime

###########
#PDF PROCESSING FOR various INPUT 
# DEPENDING ON THE PDF_PROCESSING CLASS
#by Carl Schleussner, Climate Analytics
#carl.schleussner@climateanalytics.org
###########

# IMPORT AND CONFIG 
import numpy as np
import netCDF4 as net
import dimarray as da 
import sys,glob,datetime,pickle,os,itertools
import pandas as pd
import matplotlib.pylab as plt 
from netCDF4 import Dataset,netcdftime,num2date

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
    'TXx':{'var_name':'TXX','longname':'Hot extremes (TXx)','unit':'TXx [$^\circ$ C]'},
    'TNn':{'var_name':'TNN','longname':'Cold extremes (TNn)','unit':'TNn [$^\circ$ C]'},
    'WSDI':{'var_name':'WSDI','longname':'Warm-spell duration (WSDI)','unit':'WSDI [days]'},
    'RX5':{'var_name' :'Rx5day','longname':'5-day heavy rainfall (Rx5day)','unit':'RX5 [$\%$]'},
    'RX1':{'var_name':'Rx1day','longname':'Daily heavy rainfall (Rx1day)', 'unit':'RX1 [$\%$]'}}

# time informations and periods
timeaxis=np.arange(1958,2011)
ref_period=[1960,1979]
target_periods=[[1991,2010],ref_period]
period_names=['Recent','ref']

# Set range for years for bootstrap sampling 
bs_range=[1958,2010]

# Input datasets
datasets=['HadEX2','GHCNDEX']

varoutdict={
    datasets[0]:{},    
    datasets[1]:{},    
}

# Set plottint colours
colordict={
    datasets[0]:"#247a9c",
    datasets[1]:"#df1a20",
}



read_in_data=da.read_nc('/Users/peterpfleiderer/Box Sync/0p5_observational_record/data/data_climdex/HadEx2/H2_'+varin_dict['TXx']['var_name']+'_1901-2010_RegularGrid_global_3.75x2.5deg_LSmask.nc')['Ann']
input_data=da.DimArray(read_in_data[19580101:20100101,:,:], axes=[timeaxis, read_in_data.lat, read_in_data.lon],dims=['year', 'lat', 'lon'] )

landmask=input_data.ix[10,:,:].copy()
landmask[:,:]=1
GRL_mask=Dataset('support/GRL_73x96_lat_weighted.nc4').variables['GRL'][:,:]
landmask[np.isfinite(GRL_mask)]=0 

varoutdict['GHCNDEX']['TXx']=pdf.PDF_Processing('TXx')
varoutdict['GHCNDEX']['TXx'].mask_for_ref_period_data_coverage(input_data,ref_period,check_ref_period_only=False,target_periods=target_periods,landmask=landmask,required_coverage=0.8,dataset='GHCNDEX')

varoutdict['GHCNDEX']['TXx'].derive_regional_masking(shift_lon=-180.0,region_polygons=region_polygons,mask_file='support/GHCNDEX_continent_masks.pkl')


mask=_masks['Oceania'].copy()
mask[np.isfinite(mask)]=1
varoutdict['GHCNDEX']['TXx'].plot_map(mask,
                                    color_bar=False,
                                    show=True)



