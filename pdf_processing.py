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

    def mask_for_ref_period_data_coverage(self,input_data,ref_period,maskname='global',check_ref_period_only=True,target_periods=None,landmask=None,required_coverage=None,dataset=''):
        '''
        Grid cell level averaging for different periods along time axis        
        input_data:type dimarray: dimarray with annual variable named 'data' and a 'year' -axis
        ref_period:type tuple: [startyear,endyear]
        maskname:type list : generates mask for start of reference period. Default: 'global' for first year of reference period        
        check_ref_period_only: For observations, data gaps may exist also for the reference period. False: All years investigated will be checked. 
        target_periods: Target period. Only relevant if  check_ref_period_only== True
        landmasl: np.array on lat-lon grid with 0 for ocean: This mask can be used to ignore ocean cells
        required_coverage: None or float(0-1): if None: no missings allowed. If float: ratio odf missing values per grid-cell allowed
        dataset: string: given name of the dataset
        Generates data_frame mask 
        '''    
        lat=input_data.lat
        lon=input_data.lon
        self._data=input_data

        mask_file='support/'+str(len(lat))+'x'+str(len(lon))+'_'+dataset+'_'+self._var+'_masks.pkl'

        # try to load existing mask
        if os.path.isfile(mask_file):
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

    def derive_regional_masking(self,shift_lon=0.0,regions_polygons=None,mask_file='support/144x96_srex_masks.pkl'):
        '''
        Derive regional masks
        The resulting masks can directly taken as weights for the distribution analysis. 
        shift_lon:type float: longitude shift that is required to transform the lon axis of input_data to -180 to 180
        region_polygons:type dict containing tuples: points of the polygons defining regions
        mask_file: string: location of masks (from mask_for_ref_period_data_coverage())
        '''

        input_data=self._data

        # try to load existing mask
        if os.path.isfile(mask_file):
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

        # no need to give in input_data again right?
        input_data=self._data

        # Test for time axis to be annual integers
        try: 
            timeaxis=input_data.year
        except:
            timeaxis=input_data.time 

        if isinstance(timeaxis[0], int):
            self._timeaxis=timeaxis 
        else: 
            raise ImportError("Time axis is not annual" ) 
        
        # no nead since mask is apllied later during distr
        #input_data_masked=input_data*self._masks['global']
        input_data_masked=input_data#*self._masks['global']

        # Derive time slices
        da_time_sliced=da.DimArray(axes=[np.asarray(period_names), input_data_masked.lat, input_data_masked.lon],dims=['period', 'lat', 'lon'] )        
        
        # included in the next step??
        #da_time_sliced['ref']=input_data_masked[ref_period[0]:ref_period[1]].mean(axis=0)

        for period,pname in zip(target_periods,period_names):
            print pname,period
            # why was this commented????
            #da_time_sliced[pname]=input_data_masked[period[0]:period[1]].mean(axis=0)#*self._masks[mask_for_ref_period]

            # there are nans in the period now!!!
            da_time_sliced[pname]=np.nanmean(np.array(input_data_masked[period[0]:period[1]]),axis=0)#*self._masks[mask_for_ref_period]
        
        self._data_sliced=da_time_sliced
        self._periods=self._data_sliced.period


    def derive_distributions(self,globaldist=True):
        '''
        Derive distributions for different regions. Plus regional masking.
        globaldist: type Boolean : Flag whether or not regional distributions shall be derived
        '''
        self._distributions={}
        for region in self._masks.keys():
            self._distributions[region]={}
            m=self._masks[region]
            self._distributions[region]['weight']=m[np.isfinite(m)].values.flatten()
            for key in self._periods:
                self._distributions[region][key]=self._data_sliced[key][np.isfinite(m)].values.flatten()

    def derive_pdf_difference(self,ref_period,target_period,pdf_method='python_silverman',no_hist_bins=256,range_scaling_factor=1,bin_range=None,absolute_scaling=False,relative_diff=False):

        # derive histogram pdf for pairwise differences  
        for region in self._distributions.keys():
            self._distributions[region]['pdf']={}
            self._distributions[region]['cdf']={}
            
             # get diff, relative diff is posible (in %)
            if relative_diff==False:
                diff=self._distributions[region][target_period]-self._distributions[region][ref_period]
            if relative_diff==True:
                diff=(self._distributions[region][target_period]-self._distributions[region][ref_period])/self._distributions[region][ref_period]*100

            # Set binning range for uniform analysis
            if bin_range==None:
                if absolute_scaling:
                    bin_range=[-diff.max()*range_scaling_factor,diff.max()*range_scaling_factor]
                else:
                    bin_range=[diff.min()*range_scaling_factor,diff.max()*range_scaling_factor]
            self._histogram_range=bin_range                
            

            if pdf_method=='python_silverman':
                pdf_der,no_nans=self.kernel_density_estimation(diff,bin_range,bw='silverman',kern='gaussian',region='global',method='python')
                print 'Warning, NaNs in difference kernel estimation. No of NaNs:',no_nans
                # bin x-axis is left-centered. Concert to centered 
                self._distributions[region]['pdf']['xaxis']=pdf_der[:,0]#np.asarray([(pdf_der[i,0]+pdf_der[i+1,0])/2 for i in xrange(len(pdf_der[:,0])-1)])
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
        create shuffled time slices
        '''       
        for region in self._masks.keys():
            input_data_masked=self._data.copy()[bs_range[0]:bs_range[1]]            
            mask=self._masks[region]
            self._distributions[region]['shuffled']={}
            #print input_data_masked
            for i in range(nShuff):
                # ignore nans!!
                dat = np.nanmean(input_data_masked[random.sample(input_data_masked.year,20)],axis=0).flatten()
                mdat = dat[np.where(np.isfinite(mask.flatten()))[0]]
                self._distributions[region]['shuffled'][i]=mdat


    def derive_bootstrapped_conf_interval(self,pdf_method='python_silverman',quantiles=[1,5,17,25,50,75,83,95,99],relative_diff=False):
        '''
        # derive confidence intervals for bootstrapped differences 
        quantiles: list of quantile levels (integers in [0,100])
        normed: flag, if pdf/cdf quantiles are normed to sum 1
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
                        # hist_der=np.histogram(bs_set[i]-bs_set[j],no_hist_bins, range=self._histogram_range, weights=self._distributions[region]['weight'],density=True)[0]  

                        # relative diff ???
                        if relative_diff==False:   
                            diff=bs_set[i]-bs_set[j]
                        if relative_diff==True:   
                            diff=(bs_set[i]-bs_set[j])/bs_set[j]*100

                        if pdf_method=='python_silverman':
                            hist_der,no_nans=self.kernel_density_estimation(diff,self._histogram_range,bw='silverman',kern='gaussian',region='global',method='python')
                            tot_no_nans+=no_nans
                            bs_matrix[:,index]=hist_der[:,1]#/hist_der.sum()
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


    def ks_test(self,period_name_1,period_name_2):
        '''
        2-sample KS test. Comparing distributions for 'ref' and 'recent'
        Done for each region (global is a region)
        '''
        self._ks={}
        for region in self._masks.keys():
            self._ks[region]=ks_2samp(self._distributions[region][period_name_1],self._distributions[region][period_name_2])[1]


    def kernel_density_estimation(self,diff,cutinterval,bw='silverman',kern='gaussian',region='global',method='python'):
        '''
        Fit PDFs and CDFs to respective regional functions. 
        data:dataset
        cutinterval: type tuple : Min/Max range for the PDF
        globaldist: type str : Kernel, default 'gaussian'. See R-function for update
        bw: type int : smoothing bandwidth to be used. can be a scalar or method for estimation 'scott' or 'silverman' (python only)
        region: type str : Regional analysis. 
        method: type function name : name of the function below that has to be used. default is python
        '''

        if method=='python':      
            '''
            See https://gist.github.com/tillahoffmann/f844bce2ec264c1c8cb5
            For further information on the method. 
            '''
            # passor key in self._periods:
            weights=self._distributions[region]['weight']
            # normailze weights
            weights=weights.copy()/sum(weights)

            # filter nans
            # since now np.nanmean is used for the data_slices and in bootstrapping, there shouldn't be any nans anymore
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
            # self._distributions[region][key+'_pdf']=pdf
            # cdf=pdf.copy()
            # cdf[:,1]=[sum(pdf[:i,1]) for i in xrange(pdf.shape[0])]
            # self._distributions[region][key+'_cdf']=cdf

        if method=='R':      
            '''
            Still not the same as python method
            See https://stat.ethz.ch/R-manual/R-devel/library/stats/html/density.html
            For further information on the method. 
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
                return r,np.nan
                # self._distributions[region][key+'_pdf']=r
                # cdf=r.copy()
                # cdf[:,1]=[sum(r[:i,1]) for i in xrange(r.shape[0])]
                # self._distributions[region][key+'_cdf']=cdf

    def quantiles_from_cdf(self,quantiles=[0.05,0.1,0.25,0.5,0.75,0.9,0.95]):
        '''
        estimate quantiles from cdf
        quantiles: type list : quantiles to be analyzed
        '''
        self._quantiles={}
        for region in self._masks.keys():
            self._quantiles[region]={}
            for period in self._periods:
                self._quantiles[region][period]={}
                for qu in quantiles:
                    self._quantiles[region][period][qu]=self._distributions[region][period+'_cdf'][np.argmin(abs(self._distributions[region][period+'_cdf'][:,1]-qu)),0]


    def plot_map(self,to_plot,color_bar=True,color_label=None,color_palette=plt.cm.plasma,color_range=None,limits=[-180,180,-90,90],ax=None,figsize=(8,4),coastline_width=0.5,out_file=None,title='',show=True):
        '''
        plot maps of data. 
        meta_data: list of strs: meta information required to acces data
        source: str: default='_data'. if masks are to be plotted, specify source='_masks'
        period: str: if  the averag over a period is to be plotted specify the period name
        time: int: index in time axis of data (to be plotted)
        color_bar: logical: if True, color-scale is plotted besides the map
        color_label: str: label of the color-scale
        color_palette: plt.cm. object: colors used
        color_range: [float,float]: minimal and maximal value on color-scale
        limits: [lon_min,lon_max,lat_min,lat_max]: extend of the map
        ax: subplot: subplot on which the map will be plotted
        out_file: str: location where the plot is saved
        title: str: title of the plot
        show: logical: show the subplot?
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
        #m.drawstates()
        #m.drawcountries()
        m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5) 
        m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)

        # get color_range
        if color_range==None:
            color_range=[np.min(to_plot[np.isfinite(to_plot)]),np.max(to_plot[np.isfinite(to_plot)])]

        # when saving files as pdf interpolation='none' is not working!
        #im = m.imshow(to_plot,cmap=color_palette,vmin=color_range[0],vmax=color_range[1],interpolation='nearest')

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



    def show_result(self,small_plot_function,plot_settings,output_name):
        '''
        Creates world map with little subplots around it
        Information about plot arrangement is stored in plot_settings
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
                x,y=Polygon(plot_settings['points'][region]).exterior.xy
                m.plot(x,y,'g')
                # add a point in the center of the region and a line pointing to the outersubplot
                m.plot(np.mean(x),np.mean(y),'og')
                ax_big.plot(plot_settings['line_to_subplot'][region][0:2],plot_settings['line_to_subplot'][region][2:4],'g')
                # add the outer subplot
                ax = fig.add_axes(plot_settings['subplot_window'][region],axisbg='w') 
                # fill the outer subplot with whatever is defined in small_plot_function
                small_plot_function(subax=ax,region=region)

        # add global subplot
        ax = fig.add_axes(plot_settings['subplot_window']['global'],axisbg='w') 
        small_plot_function(subax=ax,region='global')    
        
        ax_big.set_xlim([-200,200])
        ax_big.set_ylim([-100,100])

        if output_name!=None:   plt.savefig(output_name)
        if output_name==None:   plt.show()

    # def save_output(self,fname):





















