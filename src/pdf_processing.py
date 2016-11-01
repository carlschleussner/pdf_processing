"""
CLASS TO DERIVE REGIONAL AGGREGATION, PDF GENERATION AND KS TEST
FOLLOWING THE METHODOLOGY DEPLOYED
Schleussner, C.-F. et al.  ESD (2016)
by Carl Schleussner, Climate Analytics
carl.schleussner@climateanalytics.org
"""
import numpy as np
import netCDF4 as net
import sys 
sys.path.append('/home/carls/git_repos')
import dimarray as da
import itertools
import glob
import datetime
import pickle
import os 

class PDF_Processing(object):
    def __init__(self,variable_name,working_dir='./'):
        '''
        Initialize an instance of the PDF_Processing                
        variable_name:type str: Variable name 
        working_dir:type str: path to current working directory, where data would be stored as applicable 
        '''
        self._var = variable_name
        self._working_dir = working_dir

    def mask_for_ref_period_data_coverage(self,input_data,ref_period,maskname='refmask',check_ref_period_only=True,target_periods=None):
        '''
        Grid cell level averaging for different periods along time axis        
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        ref_period:type tuple: [startyear,endyear]
        maskname:type list : generates mask for start of reference period. Default: 'refmask' for first year of reference period        
        check_ref_period_only: For observations, data gaps may exist also for the reference period. False: All years investigated will be checked. 
        target_periods: Target period. Only relevant if  check_ref_period_only== True
        Generates data_frame mask 
        '''    
        if  check_ref_period_only:
            mask=np.isfinite(input_data[ref_period[0]])
            print 'No of non-NAN grid cells in Reference Mask: ', np.sum(mask)
            maskout=input_data[ref_period[0]].copy()
            maskout[:,:]=np.NaN
            maskout[mask]=1
            self._masks={maskname:maskout}

        else:
            mask=np.isfinite(input_data[ref_period[0]:ref_period[1]].mean(axis=0))
            print 'No of non-NAN grid cells in Mask over Ref period: ', np.sum(mask)            
            for tp in target_periods:
                mask*=np.isfinite(input_data[tp[0]:tp[1]].mean(axis=0))
                print 'No of non-NAN grid cells in Mask over Ref period and target period ',tp,' : ', np.sum(mask)

            maskout=input_data[ref_period[0]].copy()
            maskout[:,:]=np.NaN
            maskout[mask]=1
            self._masks={maskname:maskout}

    def derive_time_slices(self,input_data,ref_period,target_periods,period_names,mask_for_ref_period='refmask'):
        '''
        Grid cell level averaging for different periods along time axis        
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        ref_period:type tuple: [startyear,endyear]
        target_periods:type list of tuples: [[startyear,endyear],...]
        period_names:type list : names of reference periods 
        mask_for_ref_period: type str : name of mask to be deployed, default='refmask'. If set to None, no masking will be used
        Generates data_frame for each ref and target period 
        '''
        # Test for time axis to be annual integers
        timeaxis=input_data.year        
        if isinstance(input_data.year[0], int):
            self._timeaxis=timeaxis 
        else: 
            raise ImportError("Time axis is not annual" ) 
        
        input_data_masked=input_data*self._masks['refmask']
        # Derive time slices
        # period_names.append('ref')
        da_time_sliced=da.DimArray(axes=[np.asarray(period_names), input_data_masked.lat, input_data_masked.lon],dims=['period', 'lat', 'lon'] )        
        da_time_sliced['ref']=input_data_masked[ref_period[0]:ref_period[1]].mean(axis=0)

        for period,pname in zip(target_periods,period_names):
            da_time_sliced[pname]=input_data_masked[period[0]:period[1]].mean(axis=0)
        
        self._data=da_time_sliced
        self._periods=self._data.period


    def derive_regional_masking(self,regions=None):
        '''
        DUMMY FUNCTION TO BE FILLED WITH SREX REGIONAL INFORMATION. 
        '''



    def derive_distributions(self,globaldist=True,regions=None):
        '''
        Derive distributions for different regions. Plus regional masking. NOT YET IMPLEMENTED
        globaldist: type Boolean : Flag whether or not regional distributions shall be derived
        '''
        self._distributions={'global':{}}
        # GLOBAL ONLY CURRENTLY IMPLEMENTED
        for key in self._periods:            
            t=self._data[key].values.flatten()
            self._distributions['global'][key]=t[np.isfinite(t)]
        
        # Derive distribution weight (for kernel estimation)
        lat_weight=self._data['ref'].copy()
        for l in lat_weight.lon:
            lat_weight[:,l]=np.cos(np.radians(lat_weight.lat))
        exweight=lat_weight.values.flatten()
        nanfilter=np.isfinite(self._data['ref'].values.flatten())
        self._distributions['weight']=exweight[nanfilter]

    def kernel_in_R(self,cutinterval,bw,kern='gaussian',region='global'):      
        '''
        Fit PDFs and CDFs to respective regional functions. See 
        https://stat.ethz.ch/R-manual/R-devel/library/stats/html/density.html
        For further information on the method. 

        cutinterval: type tuple : Min/Max range for the PDF
        globaldist: type str : Kernel, default 'gaussian'. See R-function for update
        bw: type int : smoothing bandwidth to be used.
        region: type str : Regional analysis. Currently global only. 
        '''
        if not os.path.exists(self._working_dir+'tmp/') :
            os.mkdir(self._working_dir+'tmp/')

        for key in self._periods:
            fstring=self._working_dir+'tmp/tmp_p_to_R.dat'
            # print self._distributions[region][key].shape, self._distributions['weight'].shape,key
            out_to_r=np.vstack((self._distributions[region][key],self._distributions['weight']))
            np.savetxt(fstring,out_to_r.transpose())        
            rsavestring=self._working_dir+'tmp/tmp_processed_R_to_p.dat'

            f = open(self._working_dir+'tmp/R_kernel_ana.R', 'w') 
            f.write('t<-read.table(\"'+fstring+'\") \n')
            f.write('pdf=density(t$V1,weights=t$V2/sum(t$V2),from='+str(cutinterval[0])+',to='+str(cutinterval[1])+', kernel=\"'+kern+'\",bw='+str(bw)+',na.rm=TRUE) \n')
            f.write('out <- data.frame(x=pdf$x,y=pdf$y/sum(pdf$y))\n')
            f.write('write.table(out,\"'+rsavestring+'\" ,row.names = FALSE,col.names = FALSE ) \n')
            f.close()
            ret=os.system('Rscript '+self._working_dir+'tmp/R_kernel_ana.R')#, shell=True)
            if ret !=0:
                print 'Error in Kernel Estimation', target, model
            r=np.loadtxt(rsavestring)
            self._distributions[region][key+'_pdf']=r
            cdf=r.copy()
            cdf[:,1]=[sum(r[:i,1]) for i in xrange(r.shape[0])]
            self._distributions[region][key+'_cdf']=cdf
    
