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
        self._masks={}

    def mask_for_ref_period_data_coverage(self,input_data,ref_period,maskname='global',check_ref_period_only=True,target_periods=None,landmask=None,required_coverage=None,dataset='',overwrite=False):
        '''
        Grid cell level averaging for different periods along time axis        
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        ref_period:type tuple: [startyear,endyear]
        maskname:type list : generates mask for start of reference period. Default: 'global' for first year of reference period        
        check_ref_period_only: For observations, data gaps may exist also for the reference period. False: All years investigated will be checked. 
        target_periods: Target period. Only relevant if  check_ref_period_only== True
        landmask:type np.array on lat-lon grid with 0 for ocean: This mask can be used to ignore ocean cells
        required_coverage:type None or float(0-1): if None: no missings allowed. If float: ratio off missing values per grid-cell allowed
        dataset:type string: given name of the dataset
        Generates data_frame mask 
        '''    
        lat=input_data.lat
        lon=input_data.lon
        self._data=input_data

        mask_file='support/'+str(len(lat))+'x'+str(len(lon))+'_'+dataset+'_'+self._var+'_masks.pkl'

        # try to load existing mask
        if os.path.isfile(mask_file) and overwrite==False:
            pkl_file = open(mask_file, 'rb')
            self._masks = pickle.load(pkl_file)
            pkl_file.close()   

        # compute mask if no mask file found
        else:
            # full coverage
            # no missing values allowed
            if required_coverage==None:
                if  check_ref_period_only:
                    mask=np.isfinite(input_data[ref_period[0]:ref_period[1]].mean(axis=0))
                    print 'No of non-NAN grid cells in Reference Mask: ', np.sum(mask)

                else:
                    mask=np.isfinite(input_data[ref_period[0]:ref_period[1]].mean(axis=0))
                    print 'No of non-NAN grid cells in Mask over Ref period: ', np.sum(mask)            
                    for tp in target_periods:
                        mask*=np.isfinite(input_data[tp[0]:tp[1]].mean(axis=0))
                        print 'No of non-NAN grid cells in Mask over Ref period and target period ',tp,' : ', np.sum(mask)

            # percentage coverage
            # ratio of not missing values per grid-cell has to be larger than required_coverage
            if required_coverage!=None:
                mask=input_data[ref_period[0]].copy()
                mask[:,:]=1
                if  check_ref_period_only:
                    for y in lat:   
                        for x in lon:
                            if len(np.where(np.isfinite(input_data[ref_period[0]:ref_period[1],y,x]))[0])<input_data[ref_period[0]:ref_period[1],y,x].shape[0]*required_coverage:
                                mask[y,x]=0
                    print 'No of non-NAN grid cells in Reference Mask: ', np.sum(mask)

                else:
                    for y in lat:
                        for x in lon:
                            if len(np.where(np.isfinite(input_data[ref_period[0]:ref_period[1],y,x]))[0])<input_data[ref_period[0]:ref_period[1],y,x].shape[0]*required_coverage:
                                mask[y,x]=0
                    print 'No of non-NAN grid cells in Mask over Ref period: ', np.sum(mask)            
                    for tp in target_periods:
                        for y in lat:
                            for x in lon:
                                if len(np.where(np.isfinite(input_data[tp[0]:tp[1],y,x]))[0])<input_data[tp[0]:tp[1],y,x].shape[0]*required_coverage:
                                    mask[y,x]=0
                        print 'No of non-NAN grid cells in Mask over Ref period and target period ',tp,' : ', np.sum(mask)

            mask[mask==0]=np.nan
            mask=np.isfinite(mask)      

            # mask ocean cells
            try:
                landmask[landmask!=0]=1
                landmask[landmask==0]=np.nan
                mask*=np.isfinite(landmask)
            except:
                print 'no landmask used'


            maskout=input_data[ref_period[0]].copy()*np.nan

            # Derive land area weighting (For kernel estimation)
            lat_weight=input_data[ref_period[0]].copy()
            for l in lat_weight.lon:
                lat_weight[:,l]=np.cos(np.radians(lat_weight.lat))

            maskout[mask]=lat_weight[mask]
            # normalize 
            self._masks[maskname]=maskout/float(maskout[mask].sum())

            # save masks
            mask_output = open(mask_file, 'wb')
            pickle.dump(self._masks, mask_output)
            mask_output.close()

    def derive_regional_masking(self,shift_lon=0.0,region_polygons=None,region_type='continental',dataset='',overwrite=False):
        '''
        Derive regional masks
        The resulting masks can directly taken as weights for the distribution analysis. 
        shift_lon:type float: longitude shift that is required to transform the lon axis of input_data to -180 to 180

        # changed region_polygons!!!! for SREX i used different region_polygons. now they have to be directly Polygons!!!
        region_polygons:type dict: the polygons defining regions
        mask_file:type string: location of masks (from mask_for_ref_period_data_coverage())
        '''

        input_data=self._data
        lat=input_data.lat
        lon=input_data.lon

        mask_file='support/'+str(len(lat))+'x'+str(len(lon))+'_'+dataset+'_'+self._var+'_'+region_type+'_masks.pkl'

        # try to load existing mask
        if os.path.isfile(mask_file) and overwrite==False:
            pkl_file = open(mask_file, 'rb')
            self._masks = pickle.load(pkl_file)
            pkl_file.close()    

        else:
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
                    grid_polygons[i,j] = Polygon([(x1,y1),(x1,y2),(x2,y2),(x2,y1)])
                    #grid_polygons[i,j] = Polygon([(y1,x1),(y1,x2),(y2,x2),(y2,x1)])

            # since the lon axis has been shifted, masks and outputs will have to be shifted as well. This shift is computed here
            lon=lon-shift_lon
            shift = len(lon)-np.where(lon==lon[0]-shift_lon)[0][0]


            for region in region_polygons.keys():
                print region
                poly=region_polygons[region]
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
                    # only save mask if more than 30 grid cells available 
                    if len(np.where(np.isfinite(output))[0])>30:
                        # shift back to original longitudes
                        output=np.roll(output,shift,axis=1)
                        maskout=input_data.ix[0,:,:].copy()*np.nan
                        maskout.ix[:,:]=output[:,:]
                        self._masks[region]=maskout




                # check step. normaly one tab back
                # save masks
                mask_output = open(mask_file, 'wb')
                pickle.dump(self._masks, mask_output)
                mask_output.close()


    def derive_time_slices(self,ref_period,target_periods,period_names,mask_for_ref_period='global'):
        '''
        Grid cell level averaging for different periods along time axis        
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        ref_period:type tuple: [startyear,endyear]
        target_periods:type list of tuples: [[startyear,endyear],...]
        period_names:type list : names of reference periods 
        mask_for_ref_period: type str : name of mask to be deployed, default='global'. If set to None, no masking will be used
        Generates data_frame for each ref and target period 
        '''

        input_data=self._data

        # Load time axis (try and error). 
        try: 
            timeaxis=input_data.year
        except:
            timeaxis=input_data.time 
        if isinstance(timeaxis[0], int):
            self._timeaxis=timeaxis 
        else: 
            raise ImportError("Time axis is not annual." ) 
        

        # Derive time slices
        da_time_sliced=da.DimArray(axes=[np.asarray(period_names), input_data.lat, input_data.lon],dims=['period', 'lat', 'lon'] )        
        
        for period,pname in zip(target_periods,period_names):
            print pname,period
            da_time_sliced[pname]=np.nanmean(np.array(input_data[period[0]:period[1]]),axis=0)#*self._masks[mask_for_ref_period]
        
        self._data_sliced=da_time_sliced
        self._periods=da_time_sliced.period


    def derive_distributions(self):
        '''
        Derive distributions for different regions. Here the masks computed in mask_for_ref_period_data_coverage() and derive_regional_masking() are used.
        '''
        self._distributions={}
        for region in self._masks.keys():
            self._distributions[region]={}
            mask=self._masks[region]
            self._distributions[region]['weight']=mask[np.isfinite(mask)].values.flatten()
            for key in self._periods:
                self._distributions[region][key]=self._data_sliced[key][np.isfinite(mask)].values.flatten()

    def derive_pdf_difference(self,ref_period,target_period,pdf_method='python_silverman',bin_range=None,no_hist_bins=256,range_scaling_factor=1,absolute_scaling=False,relative_diff=False):
        '''
        Derive regional pdf's of differences between the chosen periods.
        ref_period:type str: name of the reference period
        target_period:type str: name of the target period
        pdf_method:type str: method used to derive pdf's. 'hist' for histogram; 'python_silverman' for kernel density estimation with Silverman's rule of thumb.
        bin_range:type array: bins used for kernel density estimation. If None this will be estimated in the function
        no_hist_bins:type int: number of bins used in histogram method
        absolute_scaling:type Boolean: If True, the bin_range is going to be symmetric
        range_scaling_factor:type float: factor setting the range for kernel density estimation
        relative_diff:type Boolean: if True, relative differences are considered
        '''


        # derive histogram pdf for pairwise differences  
        for region in self._distributions.keys():
            self._distributions[region]['pdf']={}
            self._distributions[region]['cdf']={}
            
             # get diff, relative diff is posible (in %)
            if relative_diff==False:
                diff=self._distributions[region][target_period]-self._distributions[region][ref_period]
            if relative_diff==True:
                diff=(self._distributions[region][target_period]-self._distributions[region][ref_period])/self._distributions[region][ref_period]*100

            # Get or set binning range for uniform analysis
            if bin_range==None:
                if absolute_scaling:
                    bin_range=[-diff.max()*range_scaling_factor,diff.max()*range_scaling_factor]
                else:
                    bin_range=[diff.min()*range_scaling_factor,diff.max()*range_scaling_factor]
            self._histogram_range=bin_range                
            

            if pdf_method=='python_silverman':
                pdf_der,no_nans=self.kernel_density_estimation(diff,bin_range,bw='silverman',region='global',method='python')
                print 'Warning, NaNs in difference kernel estimation. No of NaNs:',no_nans
                # bin x-axis is left-centered. Concert to centered 
                self._distributions[region]['pdf']['xaxis']=pdf_der[:,0]
                self._distributions[region]['cdf']['xaxis']=self._distributions[region]['pdf']['xaxis']          
                # save hist_values
                self._distributions[region]['pdf'][target_period+'_'+ref_period]=pdf_der[:,1]
                self._distributions[region]['cdf'][target_period+'_'+ref_period]=np.asarray([pdf_der[:i,1].sum() for i in xrange(len(pdf_der[:,1]))])
        
            elif pdf_method=='hist':
                der_hist=np.histogram(diff,no_hist_bins, range=self._histogram_range, weights=self._distributions[region]['weight'],density=True)
                # bin x-axis is left-centered. Concert to centered 
                self._distributions[region]['pdf']['xaxis']=np.asarray([(der_hist[1][i]+der_hist[1][i+1])/2 for i in xrange(len(der_hist[1])-1)])
                self._distributions[region]['cdf']['xaxis']=self._distributions[region]['pdf']['xaxis']          
                # save hist_values
                normed_hist=der_hist[0]/der_hist[0].sum()
                self._distributions[region]['pdf'][target_period+'_'+ref_period]=normed_hist
                self._distributions[region]['cdf'][target_period+'_'+ref_period]=np.asarray([normed_hist[:i].sum() for i in xrange(len(normed_hist))])

            else:
                print 'Method not implemented', pdf_method
                break
                

    def bootstrapping(self,bs_range,nShuff):
        '''
        Create randomly shuffled time slices
        bs_range:type list(int,int): limiting years for shuffling.
        nShuff:type int: number of required shuffled time slices
        '''       
        for region in self._masks.keys():
            input_data=self._data.copy()[bs_range[0]:bs_range[1]]            
            mask=self._masks[region]
            self._distributions[region]['shuffled']={}
            for i in range(nShuff):
                # randomly pick 20 years
                dat = np.nanmean(input_data[random.sample(input_data.year,20)],axis=0).flatten()
                mdat = dat[np.where(np.isfinite(mask.flatten()))[0]]
                self._distributions[region]['shuffled'][i]=mdat


    def derive_bootstrapped_conf_interval(self,pdf_method='python_silverman',quantiles=[1,5,17,25,50,75,83,95,99],relative_diff=False):
        '''
        Derive confidence intervals for bootstrapped differences 
        pdf_method:type str: method used to derive pdf's. 'hist' for histogram; 'python_silverman' for kernel density estimation with Silverman's rule of thumb.
        quantiles:type list of quantile levels (integers in [0,100])
        relative_diff:type Boolean: if True, relative differences are considered
        '''
        for region in self._distributions.keys():
            if self._distributions[region].has_key('pdf')==False:
                print 'Xrange not defined, run derive_pdf_difference first'
                break 
            else:
                self._distributions[region]['pdf']['bs_quantiles']={}
                self._distributions[region]['cdf']['bs_quantiles']={}
                bs_set= self._distributions[region]['shuffled']
                bs_length=len(bs_set.keys())
                no_hist_bins=len(self._distributions[region]['pdf']['xaxis'])
                bs_matrix=np.zeros((no_hist_bins,bs_length*(bs_length-1)))

                index=0
                tot_no_nans=0
                for i,j in itertools.product(xrange(bs_length),xrange(bs_length)):
                    if i != j:                           
                        # relative diff
                        if relative_diff==False:   
                            diff=bs_set[i]-bs_set[j]
                        if relative_diff==True:   
                            diff=(bs_set[i]-bs_set[j])/bs_set[j]*100

                        if pdf_method=='python_silverman':
                            hist_der,no_nans=self.kernel_density_estimation(diff,self._histogram_range,bw='silverman',region='global',method='python')
                            tot_no_nans+=no_nans
                            bs_matrix[:,index]=hist_der[:,1]
                        elif pdf_method=='hist':
                            hist_der=np.histogram(bs_set[i]-bs_set[j],no_hist_bins, range=self._histogram_range, weights=self._distributions[region]['weight'],density=True)[0]           
                            bs_matrix[:,index]=hist_der/hist_der.sum()

                        index+=1

                print 'Warning, total number of NaNs in bootstrap kernel estimation. No of NaNs:',tot_no_nans

                for qu in quantiles:
                    quant_unnormed=np.percentile(bs_matrix,qu,axis=1)
                    quant=quant_unnormed

                    self._distributions[region]['pdf']['bs_quantiles'][qu]=quant
                    self._distributions[region]['cdf']['bs_quantiles'][qu]=np.asarray([quant[:i].sum() for i in xrange(len(quant))])


    def kernel_density_estimation(self,diff,cutinterval,bw='silverman',region='global',method='python'):
        '''
        Fit PDFs and CDFs to respective regional functions. 
        data:dataset
        diff:type array: array on which the kde is applied 
        cutinterval:type tuple : Min/Max range for the PDF
        bw:type int : smoothing bandwidth to be used. can be a scalar or method for estimation 'scott' or 'silverman' (python only)
        region:type str : Regional analysis. 
        method:type function name : name of the function below that has to be used. default is python
        '''

        if method=='python':      
            '''
            See https://gist.github.com/tillahoffmann/f844bce2ec264c1c8cb5
            For further information on the method. 
            '''
            weights=self._distributions[region]['weight']
            # normailze weights
            weights=weights.copy()/sum(weights)

            # filter for nans (normally diff is coming without nans already)
            diff_nan_filter=np.isfinite(diff)
            no_nans=np.isfinite(diff).sum()-len(diff)

            # cutoff outliers
            inside_cutinterval=np.where((diff[diff_nan_filter]>=cutinterval[0]) & (diff[diff_nan_filter]<=cutinterval[1]))[0]

            kde=gaussian_kde(diff[diff_nan_filter][inside_cutinterval],weights=weights[diff_nan_filter][inside_cutinterval])
            kde.set_bandwidth(bw_method=bw)
            pdf=np.zeros([512,2])
            pdf[:,0]=np.linspace(cutinterval[0],cutinterval[1],num=512)
            pdf_zwi=kde.evaluate(pdf[:,0])
            pdf[:,1]=pdf_zwi/sum(pdf_zwi)
            return pdf,no_nans

        # comparison with a similar method in R possible:
        # if method=='R':      
        #     '''
        #     Still not the same as python method
        #     See https://stat.ethz.ch/R-manual/R-devel/library/stats/html/density.html
        #     For further information on the method. 
        #     '''
        #     if not os.path.exists(self._working_dir+'tmp/') :
        #         os.mkdir(self._working_dir+'tmp/')

        #     for key in self._periods:
        #         fstring=self._working_dir+'tmp/tmp_p_to_R.dat'
        #         out_to_r=np.vstack((self._distributions[region][key],self._distributions[region]['weight']))
        #         np.savetxt(fstring,out_to_r.transpose())        
        #         rsavestring=self._working_dir+'tmp/tmp_processed_R_to_p.dat'

        #         f = open(self._working_dir+'tmp/R_kernel_ana.R', 'w') 
        #         f.write('t<-read.table(\"'+fstring+'\") \n')
        #         f.write('pdf=density(t$V1,weights=t$V2/sum(t$V2),from='+str(cutinterval[0])+',to='+str(cutinterval[1])+', kernel=\"'+'gaussian'+'\",bw='+str(bw)+',na.rm=TRUE) \n')
        #         f.write('out <- data.frame(x=pdf$x,y=pdf$y/sum(pdf$y))\n')
        #         f.write('write.table(out,\"'+rsavestring+'\" ,row.names = FALSE,col.names = FALSE ) \n')
        #         f.close()
        #         ret=os.system('Rscript '+self._working_dir+'tmp/R_kernel_ana.R')#, shell=True)
        #         if ret !=0:
        #             print 'Error in Kernel Estimation'
        #         r=np.loadtxt(rsavestring)
        #         return r,np.nan


    def plot_map(self,to_plot,color_bar=True,color_label=None,color_palette=plt.cm.plasma,color_range=None,coastline_width=0.5,limits=[-180,180,-90,90],ax=None,figsize=(8,4),out_file=None,title='',show=True):
        '''
        Plot maps of inputted data. 
        color_bar:type logical: if True, color-scale is plotted besides the map
        color_label:type str: label of the color-scale
        color_palette:type plt.cm. object: colors used
        color_range:type [float,float]: minimal and maximal value on color-scale
        coastline_width:type float: width of coastlines
        limits:type [lon_min,lon_max,lat_min,lat_max]: extend of the map
        ax: subplot:type subplot on which the map will be plotted
        figsize:type (int,int): size of the saved figure (if saved)
        out_file:type str: location where the plot is saved. Not saved if None
        title: str:type title of the plot
        show: logical:type show the subplot?
        '''

        lat=self._data.lat.copy()
        lon=self._data.lon.copy()
        to_plot=np.array(to_plot)

        if ax==None:
            fig, ax = plt.subplots(nrows=1, ncols=1,figsize=figsize)      

        # handle 0 to 360 lon
        if max(lon)>180:
            problem_start=np.where(lon>180)[0][0]
            new_order=np.array(range(problem_start,len(lon))+range(0,problem_start))
            to_plot=to_plot[:,new_order]
            lon=lon[new_order]
            lon[lon>180]-=360


        to_plot=np.ma.masked_invalid(to_plot)

        m = Basemap(ax=ax,llcrnrlon=limits[0],urcrnrlon=limits[1],llcrnrlat=limits[2],urcrnrlat=limits[3],resolution="l",projection='cyl')
        m.drawmapboundary(fill_color='1.')

        # show coastlines and borders
        m.drawcoastlines(linewidth=coastline_width)
        m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5) 
        m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)

        # get color_range
        if color_range==None:
            color_range=[np.min(to_plot[np.isfinite(to_plot)]),np.max(to_plot[np.isfinite(to_plot)])]

        lon-=np.diff(lon,1)[0]/2.
        lat-=np.diff(lat,1)[0]/2.
        lon,lat=np.meshgrid(lon,lat)
        im = m.pcolormesh(lon,lat,to_plot,cmap=color_palette,vmin=color_range[0],vmax=color_range[1])

        # add colorbar
        if color_bar==True:
            cb = m.colorbar(im,'right', size="5%", pad="2%")
            cb.set_label(color_label, rotation=90)

        ax.set_title(title)
        ax.legend(loc='best')
        
        if out_file==None and show==True:plt.show()
        if out_file!=None:plt.tight_layout() ; plt.savefig(out_file) ; plt.clf()
        return(im)
















