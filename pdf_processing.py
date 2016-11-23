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
#sys.path.append('/home/carls/git_repos')
import dimarray as da
import itertools
import glob
import datetime
import pickle
import os 
import random as random
from weighted_kde import gaussian_kde
from six import string_types
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pylab as plt 
from scipy.stats import ks_2samp #KS-Test


class PDF_Processing(object):
    def __init__(self,variable_name,working_dir='./'):
        '''
        Initialize an instance of the PDF_Processing                
        variable_name:type str: Variable name 
        working_dir:type str: path to current working directory, where data would be stored as applicable 
        '''
        self._var = variable_name
        self._working_dir = working_dir

    def mask_for_ref_period_data_coverage(self,input_data,ref_period,maskname='global',check_ref_period_only=True,target_periods=None):
        '''
        Grid cell level averaging for different periods along time axis        
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        ref_period:type tuple: [startyear,endyear]
        maskname:type list : generates mask for start of reference period. Default: 'global' for first year of reference period        
        check_ref_period_only: For observations, data gaps may exist also for the reference period. False: All years investigated will be checked. 
        target_periods: Target period. Only relevant if  check_ref_period_only== True
        Generates data_frame mask 
        '''    
        if  check_ref_period_only:
            mask=np.isfinite(input_data[ref_period[0]])
            print 'No of non-NAN grid cells in Reference Mask: ', np.sum(mask)

        else:
            mask=np.isfinite(input_data[ref_period[0]:ref_period[1]].mean(axis=0))
            print 'No of non-NAN grid cells in Mask over Ref period: ', np.sum(mask)            
            for tp in target_periods:
                mask*=np.isfinite(input_data[tp[0]:tp[1]].mean(axis=0))
                print 'No of non-NAN grid cells in Mask over Ref period and target period ',tp,' : ', np.sum(mask)

        maskout=input_data[ref_period[0]].copy()
        maskout[:,:]=np.NaN

        # Derive distribution weight (For kernel estimation)
        lat_weight=input_data[ref_period[0]].copy()
        for l in lat_weight.lon:
            lat_weight[:,l]=np.cos(np.radians(lat_weight.lat))

        maskout[mask]=lat_weight[mask]
        self._masks={maskname:np.ma.masked_invalid(maskout)}



    def derive_regional_masking(self,input_data,shift_lon=0.0,regions_polygons=None):
        '''
        Derive regional masks
        The resulting masks can directly taken as weights for the distribution analysis. 
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        shift_lon:type float: longitude shift that is required to transform the lon axis of input_data to -180 to 180
        region_polygons:type dict containing tuples: points of the polygons defining regions
        '''

        # get information about grid of input data
        lat = input_data.lat
        lon = input_data.lon.squeeze()+shift_lon
        nx = len(lon)
        ny = len(lat)

        # loop over the grid to get grid polygons
        grid_polygons = np.empty((nx,ny),dtype=Polygon)
        dx = np.zeros((nx))
        dy = np.zeros((ny))
        dx[1:] = np.abs(np.diff(lon,1))
        dx[0] = dx[1]
        dy[1:] = np.abs(np.diff(lat,1))
        dy[0] = dy[1]
        for i in range(nx):
            x1 = lon[i]-dx[i]/2.
            x2 = lon[i]+dx[i]/2.
            for j in range(ny):
                y1 = lat[j]-dy[j]/2.
                y2 = lat[j]+dy[j]/2.
                #grid_polygons[i,j] = Polygon([(x1,y1),(x1,y2),(x2,y2),(x2,y1)])
                grid_polygons[i,j] = Polygon([(y1,x1),(y1,x2),(y2,x2),(y2,x1)])

        # since the lon axis has been shifted, masks and outputs will have to be shifted as well. This shift is computed here
        lon=lon-shift_lon
        shift = len(lon)-np.where(lon==lon[0]-shift_lon)[0][0]


        for region in regions_polygons.keys():
            poly=Polygon(regions_polygons[region]['poly'])
            overlap = np.zeros((ny,nx))
            for i in range(nx):
                for j in range(ny):
                    # check whether data exists in grid cell
                    if np.isfinite(np.roll(self._masks['global'],shift,axis=1)[j,i]):
                        # get fraction of grid-cell covered by polygon
                        intersect = grid_polygons[i,j].intersection(poly).area/grid_polygons[i,j].area*poly.area
                        # multiply overlap with latitude weighting
                        overlap[j,i] = intersect*np.cos(np.radians(lat[j]))

            # renormalize overlap to get sum(mask)=1
            overlap_zwi=overlap.copy()
            overlap_sum=sum(overlap_zwi.flatten())
            if overlap_sum!=0:
                output=np.zeros(overlap.shape)
                output=overlap/overlap_sum
                # mask zeros
                output[output==0]=np.nan
                output=np.ma.masked_invalid(output)
                # only save mask if more than 50 grid cells available 
                if len(np.where(np.isfinite(output))[0])>50:
                    # shift back to original longitudes
                    self._masks[region]=np.roll(output,shift,axis=1)

    def derive_time_slices(self,input_data,ref_period,target_periods,period_names,mask_for_ref_period='global'):
        '''
        Grid cell level averaging for different periods along time axis        
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        ref_period:type tuple: [startyear,endyear]
        target_periods:type list of tuples: [[startyear,endyear],...]
        period_names:type list : names of reference periods 
        mask_for_ref_period: type str : name of mask to be deployed, default='global'. If set to None, no masking will be used
        Generates data_frame for each ref and target period 
        '''
        # Test for time axis to be annual integers
        timeaxis=input_data.year        
        if isinstance(input_data.year[0], int):
            self._timeaxis=timeaxis 
        else: 
            raise ImportError("Time axis is not annual" ) 
        
        # no nead since mask is apllied later during distr
        #input_data_masked=input_data*self._masks['global']
        input_data_masked=input_data#*self._masks['global']

        # Derive time slices
        # period_names.append('ref')
        da_time_sliced=da.DimArray(axes=[np.asarray(period_names), input_data_masked.lat, input_data_masked.lon],dims=['period', 'lat', 'lon'] )        
        da_time_sliced['ref']=input_data_masked[ref_period[0]:ref_period[1]].mean(axis=0)

        for period,pname in zip(target_periods,period_names):
            print pname,period
            da_time_sliced[pname]=input_data_masked[period[0]:period[1]].mean(axis=0)
        
        self._data=da_time_sliced
        self._periods=self._data.period


    def derive_distributions(self,globaldist=True):
        '''
        Derive distributions for different regions. Plus regional masking. NOT YET IMPLEMENTED
        globaldist: type Boolean : Flag whether or not regional distributions shall be derived
        '''
        self._distributions={}

        # REGIONAL (SREX regions)
        # mask grid-cells where no data exists in reference
        # do we still need this nan filter??? isnt mask there for that?
        nanfilter=np.isfinite(self._data['ref'].values.flatten())
        for region in self._masks.keys():
            self._distributions[region]={}
            # use mask calculated in derive_regional_masking()
            mask=np.isfinite(self._masks[region].flatten()[nanfilter])
            self._distributions[region]['weight']=self._masks[region].flatten()[nanfilter][mask] 
            for key in self._periods:
                t=self._data[key].values.flatten()
                self._distributions[region][key]=t[nanfilter][mask]   

    def bootstrapping(self,input_data,nShuff):
        # Test for time axis to be annual integers
        timeaxis=input_data.year        
        if isinstance(input_data.year[0], int):
            self._timeaxis=timeaxis 
        else: 
            raise ImportError("Time axis is not annual" ) 
        
        # do we still need this nan filter??? isnt mask there for that?
        nanfilter=np.isfinite(self._data['ref'].values.flatten())

        for region in self._masks.keys():
            input_data_masked=input_data
            mask=np.isfinite(self._masks[region].flatten()[nanfilter])
            self._distributions[region]['shuffled']={}
            for i in range(nShuff):
                dat = input_data_masked[random.sample(timeaxis,20)]
                mdat = np.ma.masked_array(dat,np.isnan(dat))
                # input_sample is comparable to da_time_sliced
                # .filled(np.nan) is required to maintain the same dimensions as ._distributions[region]['weight']
                input_sample=np.mean(mdat,axis=0).filled(np.nan)
                t=input_sample.flatten()
                self._distributions[region]['shuffled'][i]=t[nanfilter][mask]              

    def simple_difference(self):
        '''
        Get the difference in index between 'Recent' and 'ref'
        Bootstrapping difference + quantiles
        '''
        self._mean={}
        self._difference={}
        for region in self._distributions.keys():
            self._mean[region]={}
            for period in self._periods:
                # compute weighted mean for each period
                self._mean[region][period]=np.sum(self._distributions[region][period]*self._distributions[region]['weight'])/np.sum(self._distributions[region]['weight'])

            self._mean[region]['shuffled']=[]
            for i in self._distributions[region]['shuffled'].keys():
                # computed weighted mean for each shuffled sample
                self._mean[region]['shuffled'].append(np.sum(self._distributions[region]['shuffled'][i]*self._distributions[region]['weight'])/np.sum(self._distributions[region]['weight']))

            self._difference[region]={}
            self._difference[region]['Recent-ref']=self._mean[region]['Recent']-self._mean[region]['ref']
            self._difference[region]['quantiles']=[]
            # compute quantiles for differences between ref and shuffled realizations
            for qu in [1,5,10,25,50,75,90,95,99]:
                self._difference[region]['quantiles'].append((qu,np.percentile(self._mean[region]['shuffled']-self._mean[region]['ref'],qu)))

    def ks_test(self):
        '''
        2-sample KS test. Comparing distributions for 'ref' and 'recent'
        Done for each region (global is a region)
        '''
        self._ks={}
        for region in self._masks.keys():
            self._ks[region]=ks_2samp(self._distributions[region]['ref'],self._distributions[region]['Recent'])[1]

    def kernel_in_R(self,cutinterval,bw,kern='gaussian',region='global'):      
        '''
        Fit PDFs and CDFs to respective regional functions. See 
        https://stat.ethz.ch/R-manual/R-devel/library/stats/html/density.html
        For further information on the method. 
        cutinterval: type tuple : Min/Max range for the PDF
        globaldist: type str : Kernel, default 'gaussian'. See R-function for update
        bw: type int : smoothing bandwidth to be used.
        region: type str : Regional analysis. 
        '''
        if not os.path.exists(self._working_dir+'tmp/') :
            os.mkdir(self._working_dir+'tmp/')

        for key in self._periods:
            fstring=self._working_dir+'tmp/tmp_p_to_R.dat'
            out_to_r=np.vstack((self._distributions[region][key],self._distributions[region]['weight']))
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
                print 'Error in Kernel Estimation'
            r=np.loadtxt(rsavestring)
            self._distributions[region][key+'_pdf']=r
            cdf=r.copy()
            cdf[:,1]=[sum(r[:i,1]) for i in xrange(r.shape[0])]
            self._distributions[region][key+'_cdf']=cdf

    def kernel_in_PY(self,cutinterval,bw,kern='gaussian',region='global'):      
        '''
        Fit PDFs and CDFs to respective regional functions. See 
        https://gist.github.com/tillahoffmann/f844bce2ec264c1c8cb5
        For further information on the method. 
        cutinterval: type tuple : Min/Max range for the PDF
        globaldist: type str : Kernel, default 'gaussian'. See R-function for update
        bw: type int : smoothing bandwidth to be used. can be a scalar or method for estimation 'scott' or 'silverman'
        region: type str : Regional analysis. 
        '''
        if not os.path.exists(self._working_dir+'tmp/') :
            os.mkdir(self._working_dir+'tmp/')

        for key in self._periods:
            weights=self._distributions[region]['weight']
            # normailze weights
            weights=weights.copy()/sum(weights)
            kde=gaussian_kde(self._distributions[region][key],weights=weights)
            kde.set_bandwidth(bw_method=bw)
            pdf=np.zeros([512,2])
            pdf[:,0]=np.linspace(cutinterval[0],cutinterval[1],num=512)
            pdf_zwi=kde.evaluate(pdf[:,0])
            pdf[:,1]=pdf_zwi/sum(pdf_zwi)
            self._distributions[region][key+'_pdf']=pdf
            cdf=pdf.copy()
            cdf[:,1]=[sum(pdf[:i,1]) for i in xrange(pdf.shape[0])]
            self._distributions[region][key+'_cdf']=cdf


    def show_maps(self,output_name):
        '''
        Creates map with information on grid-point level
        '''
        fig=plt.figure(figsize=(6,4))
        m=Basemap(llcrnrlon=-180,urcrnrlon=180,llcrnrlat=-90,urcrnrlat=90,resolution="l",projection='cyl')
        m.drawcoastlines()
        Z=self._data['Recent']-self._data['ref']

        # otherwise there is a bug, don't know how to fix this (this gridcell is on the edge)
        Z[:,180]=np.nan
        Zm=np.ma.masked_invalid(Z)

        # lon and lat have to be meshgrid and they should indicate the edges of gridcells, not centers of grid cells
        lon=self._data.lon.copy()
        lon[lon>180]-=360
        lon-=3.75/2
        lon=np.append(lon,np.array(lon[-1]+3.75))

        lat=self._data.lat.copy()
        lat-=2.5/2
        lat=np.append(lat,lat[-1]+2.5)

        lons, lats = np.meshgrid(lon,lat)

        im1 = m.pcolormesh(lons,lats,Zm,cmap='PiYG',vmin=-5,vmax=5)
        fig.colorbar(im1,orientation='horizontal')

        if output_name!=None:   plt.savefig(output_name)
        if output_name==None:   plt.show()


    def show_result(self,small_plot_function,srex_polygons,output_name):
        '''
        Creates world map with little subplots around it
        Information about plot arrangement is stored in srex_polygons
        The subplots will be filled with the input of small_plot_function
        small_plot_functions can be defined outside of this class. the function name of the defined small_plot_function is given to show_results()
        '''

        # settings for big plot image
        ratio=0.2
        fig = plt.figure(figsize=(9,6))

        # big plot window
        ax_big=fig.add_axes([0,0,1,1])
        ax_big.axis('off')

        # map in the center of the big plot window (transparant background)
        ax_map=fig.add_axes([ratio,ratio,1-2*ratio,1-2*ratio])
        ax_map.patch.set_facecolor('None')
        ax_map.axis('off')
        m=Basemap(ax=ax_map)
        m.drawcoastlines()

        for region in self._masks.keys():
            if region != 'global':
                # plot the ploygon on the map
                x,y=Polygon(srex_polygons[region]['points']).exterior.xy
                m.plot(x,y,'g')
                # add a point in the center of the region and a line pointing to the outersubplot
                m.plot(np.mean(x),np.mean(y),'og')
                ax_big.plot(srex_polygons[region]['line-out_plot'][0:2],srex_polygons[region]['line-out_plot'][2:4],'g')
                # add the outer subplot
                ax = fig.add_axes(srex_polygons[region]['out_plot'],axisbg='w') 
                # fill the outer subplot with whatever is defined in small_plot_function
                small_plot_function(subax=ax,region=region)

        # add global subplot
        ax = fig.add_axes([0.01,0.005,0.21,0.28],axisbg='w') 
        small_plot_function(subax=ax,region='global')    
        
        ax_big.set_xlim([-200,200])
        ax_big.set_ylim([-100,100])

        if output_name!=None:   plt.savefig(output_name)
        if output_name==None:   plt.show()























