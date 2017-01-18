#############
# SREX plot info
#############

import numpy as np
import pickle
from shapely.geometry import Polygon, MultiPolygon


#############
# all regions
#############

pkl_file = open('srex_dict.pkl', 'rb')
srex_polygons = pickle.load(pkl_file)
pkl_file.close()    

plot_settings={'subplot_window':{}, 'line_to_subplot':{}, 'points':{}}
	
plot_settings['subplot_window']['CAM']=[0/10.,4/6.,0.09,1/6.]
plot_settings['subplot_window']['AMZ']=[0/10.,3/6.,0.09,1/6.]
plot_settings['subplot_window']['NEB']=[0/10.,2/6.,0.09,1/6.]
plot_settings['subplot_window']['WSA']=[0/10.,1/6.,0.09,1/6.]

plot_settings['subplot_window']['WNA']=[0/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['ALA']=[1/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['CNA']=[2/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['ENA']=[3/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['CGI']=[4/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['NEU']=[5/10.,5/6.,0.09,1/6.] 
plot_settings['subplot_window']['CEU']=[6/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['CAS']=[7/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['NAS']=[8/10.,5/6.,0.09,1/6.]
plot_settings['subplot_window']['TIB']=[9/10.,5/6.,0.09,1/6.]
	
plot_settings['subplot_window']['EAS']=[9/10.,4/6.,0.09,1/6.]
plot_settings['subplot_window']['SAS']=[9/10.,3/6.,0.09,1/6.]
plot_settings['subplot_window']['NAU']=[9/10.,2/6.,0.09,1/6.]
plot_settings['subplot_window']['SEA']=[9/10.,1/6.,0.09,1/6.]

# plot_settings['subplot_window']['WNA']=[0/10.,0/6.,0.09,1/6.]
# plot_settings['subplot_window']['SAU']=[1/10.,0/6.,0.09,1/6.]
plot_settings['subplot_window']['SSA']=[2/10.,0/6.,0.09,1/6.]
plot_settings['subplot_window']['MED']=[3/10.,0/6.,0.09,1/6.]
plot_settings['subplot_window']['WAF']=[4/10.,0/6.,0.09,1/6.]
plot_settings['subplot_window']['SAH']=[5/10.,0/6.,0.09,1/6.] 
plot_settings['subplot_window']['SAF']=[6/10.,0/6.,0.09,1/6.]
plot_settings['subplot_window']['EAF']=[7/10.,0/6.,0.09,1/6.]
plot_settings['subplot_window']['WAS']=[8/10.,0/6.,0.09,1/6.]
plot_settings['subplot_window']['SAU']=[9/10.,0/6.,0.09,1/6.]

x_ratio=0.67
y_ratio=0.5

for region in srex_polygons.keys():
	if region not in ['global']:
		# stor original polygon information
		plot_settings['points'][region]=srex_polygons[region]['points']

		x,y=Polygon(srex_polygons[region]['points']).exterior.xy
		x_mean,y_mean=np.mean(x),np.mean(y)
		x_out,y_out=plot_settings['subplot_window'][region][0],plot_settings['subplot_window'][region][1]

		if x_out==0.0 and y_out==0.0: x_out+=1/10.;y_out+=1/6.
		elif x_out==9/10. and y_out==5/6.: x_out+=0;y_out+=0
		elif x_out==9/10. and y_out==0.0: x_out+=0.0;y_out+=1/6.
		elif x_out==0.0 and y_out==5/6.: x_out+=1/10.;y_out+=0

		elif x_out==0.0: x_out+=1/18.;y_out+=1/20.
		elif x_out==9/10.: x_out+=0.0;y_out+=1/20.
		elif y_out==0.0: x_out+=1/18.;y_out+=1/6.
		elif y_out==5/6.: x_out+=1/18.;y_out+=0
		plot_settings['line_to_subplot'][region]=[x_mean*x_ratio,(x_out)*400-200,y_mean*y_ratio,(y_out)*200-100]

for region in plot_settings['subplot_window'].keys():
	plot_settings['subplot_window'][region][0]+=0.01
	plot_settings['subplot_window'][region][1]+=0.005
	plot_settings['subplot_window'][region][2]+=0.005
	plot_settings['subplot_window'][region][3]-=0.06

plot_settings['subplot_window']['global']=[0.01,0.005,0.095*2,1/6.-0.06]

output = open('plot_settings.pkl', 'wb')
pickle.dump(plot_settings, output)
output.close()


# #############
# # all regions
# #############

# pkl_file = open('srex_dict.pkl', 'rb')
# srex_polygons = pickle.load(pkl_file)
# pkl_file.close()    

	
# srex_polygons['CAM']['out_plot']=[0/9.,4/6.,0.1,1/6.]
# srex_polygons['AMZ']['out_plot']=[0/9.,3/6.,0.1,1/6.]
# srex_polygons['WSA']['out_plot']=[0/9.,2/6.,0.1,1/6.]
# #srex_polygons['SEA']['out_plot']=[0/9.,1/6,0.1,1/6.]

# srex_polygons['ALA']['out_plot']=[0/9.,5/6.,0.1,1/6.]
# srex_polygons['WNA']['out_plot']=[1/9.,5/6.,0.1,1/6.]
# srex_polygons['CNA']['out_plot']=[2/9.,5/6.,0.1,1/6.]
# srex_polygons['ENA']['out_plot']=[3/9.,5/6.,0.1,1/6.]
# srex_polygons['CGI']['out_plot']=[4/9.,5/6.,0.1,1/6.]
# srex_polygons['CEU']['out_plot']=[5/9.,5/6.,0.1,1/6.] 
# srex_polygons['NAS']['out_plot']=[6/9.,5/6.,0.1,1/6.]
# srex_polygons['CAS']['out_plot']=[7/9.,5/6.,0.1,1/6.]
# srex_polygons['TIB']['out_plot']=[8/9.,5/6.,0.1,1/6.]
	
# srex_polygons['EAS']['out_plot']=[8/9.,4/6.,0.1,1/6.]
# srex_polygons['SAS']['out_plot']=[8/9.,3/6.,0.1,1/6.]
# srex_polygons['NAU']['out_plot']=[8/9.,2/6.,0.1,1/6.]
# srex_polygons['SAU']['out_plot']=[8/9.,1/6.,0.1,1/6.]

# #srex_polygons['WNA']['out_plot']=[0/9.,0/6.,0.1,1/6.]
# #srex_polygons['ALA']['out_plot']=[1/9.,0/6.,0.1,1/6.]
# srex_polygons['SSA']['out_plot']=[2/9.,0/6.,0.1,1/6.]
# srex_polygons['NEB']['out_plot']=[3/9.,0/6.,0.1,1/6.]
# srex_polygons['MED']['out_plot']=[4/9.,0/6.,0.1,1/6.]
# srex_polygons['SAH']['out_plot']=[5/9.,0/6.,0.1,1/6.] 
# srex_polygons['WAF']['out_plot']=[6/9.,0/6.,0.1,1/6.]
# srex_polygons['SAF']['out_plot']=[7/9.,0/6.,0.1,1/6.]
# srex_polygons['EAF']['out_plot']=[8/9.,0/6.,0.1,1/6.]

# srex_polygons['SAU']['out_plot']=[8/9.,1/6.,0.1,1/6.]
# srex_polygons['SAU']['out_plot']=[8/9.,1/6.,0.1,1/6.]
# srex_polygons['SAU']['out_plot']=[8/9.,1/6.,0.1,1/6.]
# srex_polygons['SAU']['out_plot']=[8/9.,1/6.,0.1,1/6.]



# x_ratio=0.67
# y_ratio=0.5

# for region in srex_polygons.keys():
# 	if region not in ['global']:
# 		x,y=Polygon(srex_polygons[region]['points']).exterior.xy
# 		x_mean,y_mean=np.mean(x),np.mean(y)
# 		x_out,y_out=srex_polygons[region]['out_plot'][0],srex_polygons[region]['out_plot'][1]

# 		if x_out==0.0 and y_out==0.0: x_out+=1/9.;y_out+=1/6.
# 		elif x_out==8/9. and y_out==5/6.: x_out+=0;y_out+=0
# 		elif x_out==8/9. and y_out==0.0: x_out+=0.0;y_out+=1/6.
# 		elif x_out==0.0 and y_out==5/6.: x_out+=1/9.;y_out+=0

# 		elif x_out==0.0: x_out+=1/18.;y_out+=1/12.
# 		elif x_out==8/9.: x_out+=0.0;y_out+=1/12.
# 		elif y_out==0.0: x_out+=1/18.;y_out+=1/6.
# 		elif y_out==5/6.: x_out+=1/18.;y_out+=0
# 		srex_polygons[region]['line-out_plot']=[x_mean*x_ratio,(x_out)*400-200,y_mean*y_ratio,(y_out)*200-100]


# for region in srex_polygons.keys():
# 	srex_polygons[region]['out_plot'][0]+=0.01
# 	srex_polygons[region]['out_plot'][1]+=0.005
# 	srex_polygons[region]['out_plot'][2]+=0.005
# 	srex_polygons[region]['out_plot'][3]-=0.06

# #del srex_polygons['global']

# output = open('srex_dict_all.pkl', 'wb')
# pickle.dump(srex_polygons, output)
# output.close()

################
# reduced form
##################

# pkl_file = open('srex_dict.pkl', 'rb')
# srex_polygons = pickle.load(pkl_file)
# pkl_file.close()    

# srex_polygons['WNA']['out_plot']=[0.0,0.6,0.1,0.2]
# srex_polygons['CAM']['out_plot']=[0.0,0.4,0.1,0.2]
# srex_polygons['AMZ']['out_plot']=[0.0,0.2,0.1,0.2]

# srex_polygons['WNA']['out_plot']=[0.0,0.8,0.1,0.2]
# srex_polygons['ALA']['out_plot']=[0.125,0.8,0.1,0.2]
# srex_polygons['CNA']['out_plot']=[0.25,0.8,0.1,0.2]
# srex_polygons['CGI']['out_plot']=[0.375,0.8,0.1,0.2]
# srex_polygons['NEU']['out_plot']=[0.5,0.8,0.1,0.2]
# srex_polygons['CEU']['out_plot']=[0.625,0.8,0.1,0.2] 
# srex_polygons['NAS']['out_plot']=[0.75,0.8,0.1,0.2]
# srex_polygons['CAS']['out_plot']=[0.875,0.8,0.1,0.2]
	
# srex_polygons['TIB']['out_plot']=[0.875,0.6,0.1,0.2]
# srex_polygons['EAS']['out_plot']=[0.875,0.4,0.1,0.2]
# srex_polygons['SEA']['out_plot']=[0.875,0.2,0.1,0.2]

# srex_polygons['WSA']['out_plot']=[0.0,0.0,0.1,0.2]
# srex_polygons['SSA']['out_plot']=[0.125,0.0,0.1,0.2]
# srex_polygons['MED']['out_plot']=[0.25,0.0,0.1,0.2]
# srex_polygons['SAH']['out_plot']=[0.375,0.0,0.1,0.2]
# srex_polygons['WAS']['out_plot']=[0.5,0.0,0.1,0.2]
# srex_polygons['SAS']['out_plot']=[0.625,0.0,0.1,0.2]
# srex_polygons['SAU']['out_plot']=[0.75,0.0,0.1,0.2]
# srex_polygons['NAU']['out_plot']=[0.875,0.0,0.1,0.2]



# x_ratio=0.67
# y_ratio=0.5

# for region in srex_polygons.keys():
# 	if region not in ['EAF','SAF','WAF','global']:
# 		x,y=Polygon(srex_polygons[region]['points']).exterior.xy
# 		x_mean,y_mean=np.mean(x),np.mean(y)
# 		x_out,y_out=srex_polygons[region]['out_plot'][0],srex_polygons[region]['out_plot'][1]

# 		if x_out==0.0 and y_out==0.0: x_out+=0.125;y_out+=0.2
# 		elif x_out==0.875 and y_out==0.8: x_out+=0;y_out+=0
# 		elif x_out==0.875 and y_out==0.0: x_out+=0.0;y_out+=0.2
# 		elif x_out==0.0 and y_out==0.8: x_out+=0.125;y_out+=0

# 		elif x_out==0.0: x_out+=0.0625;y_out+=0.1
# 		elif x_out==0.875: x_out+=0.0;y_out+=0.1
# 		elif y_out==0.0: x_out+=0.0625;y_out+=0.2
# 		elif y_out==0.8: x_out+=0.0625;y_out+=0
# 		srex_polygons[region]['line-out_plot']=[x_mean*x_ratio,(x_out)*400-200,y_mean*y_ratio,(y_out)*200-100]


# for region in srex_polygons.keys():
# 	srex_polygons[region]['out_plot'][0]+=0.01
# 	srex_polygons[region]['out_plot'][1]+=0.005
# 	srex_polygons[region]['out_plot'][2]+=0.005
# 	srex_polygons[region]['out_plot'][3]-=0.06

# #del srex_polygons['global']

# output = open('srex_dict_small.pkl', 'wb')
# pickle.dump(srex_polygons, output)
# output.close()



# pkl_file = open('srex_dict.pkl', 'rb')
# srex_polygons = pickle.load(pkl_file)
# pkl_file.close()    

# srex_polygons['WNA']['out_plot']=[0.0,0.6,0.1,0.2]
# srex_polygons['CAM']['out_plot']=[0.0,0.4,0.1,0.2]
# srex_polygons['AMZ']['out_plot']=[0.0,0.2,0.1,0.2]

# srex_polygons['ALA']['out_plot']=[0.0,0.8,0.1,0.2]
# srex_polygons['CNA']['out_plot']=[0.125,0.8,0.1,0.2]
# srex_polygons['ENA']['out_plot']=[0.25,0.8,0.1,0.2]
# srex_polygons['CGI']['out_plot']=[0.375,0.8,0.1,0.2]
# srex_polygons['NEU']['out_plot']=[0.5,0.8,0.1,0.2]
# srex_polygons['CEU']['out_plot']=[0.625,0.8,0.1,0.2] 
# srex_polygons['NAS']['out_plot']=[0.75,0.8,0.1,0.2]
# srex_polygons['CAS']['out_plot']=[0.875,0.8,0.1,0.2]
	
# srex_polygons['WSA']['out_plot']=[0.0,0.0,0.1,0.2]
# srex_polygons['SSA']['out_plot']=[0.125,0.0,0.1,0.2]
# srex_polygons['MED']['out_plot']=[0.25,0.0,0.1,0.2]
# srex_polygons['SAH']['out_plot']=[0.375,0.0,0.1,0.2]
# srex_polygons['WAS']['out_plot']=[0.5,0.0,0.1,0.2]
# srex_polygons['SAS']['out_plot']=[0.625,0.0,0.1,0.2]
# srex_polygons['SAU']['out_plot']=[0.75,0.0,0.1,0.2]
# srex_polygons['NAU']['out_plot']=[0.875,0.0,0.1,0.2]

# srex_polygons['TIB']['out_plot']=[0.875,0.6,0.1,0.2]
# srex_polygons['EAS']['out_plot']=[0.875,0.4,0.1,0.2]
# srex_polygons['SEA']['out_plot']=[0.875,0.2,0.1,0.2]

# x_ratio=0.67
# y_ratio=0.5

# for region in srex_polygons.keys():
# 	if region not in ['EAF','SAF','WAF']:
# 		x,y=Polygon(srex_polygons[region]['points']).exterior.xy
# 		x_mean,y_mean=np.mean(x),np.mean(y)
# 		x_out,y_out=srex_polygons[region]['out_plot'][0],srex_polygons[region]['out_plot'][1]

# 		if x_out==0.0 and y_out==0.0: x_out+=0.125;y_out+=0.2
# 		elif x_out==0.875 and y_out==0.8: x_out+=0;y_out+=0
# 		elif x_out==0.875 and y_out==0.0: x_out+=0.0;y_out+=0.2
# 		elif x_out==0.0 and y_out==0.8: x_out+=0.125;y_out+=0

# 		elif x_out==0.0: x_out+=0.0625;y_out+=0.1
# 		elif x_out==0.875: x_out+=0.0;y_out+=0.1
# 		elif y_out==0.0: x_out+=0.0625;y_out+=0.2
# 		elif y_out==0.8: x_out+=0.0625;y_out+=0
# 		srex_polygons[region]['line-out_plot']=[x_mean*x_ratio,(x_out)*400-200,y_mean*y_ratio,(y_out)*200-100]


# for region in srex_polygons.keys():
# 	srex_polygons[region]['out_plot'][0]+=0.01
# 	srex_polygons[region]['out_plot'][1]+=0.005
# 	srex_polygons[region]['out_plot'][2]+=0.005
# 	srex_polygons[region]['out_plot'][3]-=0.06

# output = open('srex_dict.pkl', 'wb')
# pickle.dump(srex_polygons, output)
# output.close()




# pkl_file = open('srex_dict.pkl', 'rb')
# srex_polygons = pickle.load(pkl_file)
# pkl_file.close()    

# srex_polygons['WNA']['out_plot']=[0,0.6,0.1,0.2]
# srex_polygons['CAM']['out_plot']=[0,0.4,0.1,0.2]
# srex_polygons['AMZ']['out_plot']=[0,0.2,0.1,0.2]

# srex_polygons['ALA']['out_plot']=[0.1,0.8,0.1,0.2]
# srex_polygons['CNA']['out_plot']=[0.2,0.8,0.1,0.2]
# srex_polygons['ENA']['out_plot']=[0.3,0.8,0.1,0.2]
# srex_polygons['CGI']['out_plot']=[0.4,0.8,0.1,0.2]
# srex_polygons['NEU']['out_plot']=[0.5,0.8,0.1,0.2] 
# srex_polygons['CEU']['out_plot']=[0.6,0.8,0.1,0.2]
# srex_polygons['NAS']['out_plot']=[0.7,0.8,0.1,0.2]
# srex_polygons['CAS']['out_plot']=[0.8,0.8,0.1,0.2]
# srex_polygons['TIB']['out_plot']=[0.9,0.8,0.1,0.2]

# srex_polygons['WSA']['out_plot']=[0.1,0.0,0.1,0.2]
# srex_polygons['NEB']['out_plot']=[0.2,0.0,0.1,0.2]
# srex_polygons['SSA']['out_plot']=[0.3,0.0,0.1,0.2]
# srex_polygons['SAH']['out_plot']=[0.4,0.0,0.1,0.2]
# srex_polygons['MED']['out_plot']=[0.5,0.0,0.1,0.2]
# srex_polygons['WAS']['out_plot']=[0.6,0.0,0.1,0.2]
# srex_polygons['SAS']['out_plot']=[0.7,0.0,0.1,0.2]
# srex_polygons['SAU']['out_plot']=[0.8,0.0,0.1,0.2]

# srex_polygons['EAS']['out_plot']=[0.9,0.6,0.1,0.2]
# srex_polygons['SEA']['out_plot']=[0.9,0.4,0.1,0.2]
# srex_polygons['NAU']['out_plot']=[0.9,0.2,0.1,0.2]

# for region in srex_polygons.keys():
# 	if region not in ['EAF','SAF','WAF']:
# 		x,y=Polygon(srex_polygons[region]['points']).exterior.xy
# 		x_mean,y_mean=np.mean(x),np.mean(y)
# 		x_out,y_out=srex_polygons[region]['out_plot'][0],srex_polygons[region]['out_plot'][1]

# 		if x_out==0.0 and y_out==0.0:srex_polygons[region]['line-out_plot']=[x_mean*0.60,(x_out+0.1)*400-200,y_mean*0.6,(y_out+0.2)*200-100]
# 		elif x_out==0.9 and y_out==0.8:srex_polygons[region]['line-out_plot']=[x_mean*0.60,(x_out)*400-200,y_mean*0.6,(y_out)*200-100]
# 		elif x_out==0.9 and y_out==0.0:srex_polygons[region]['line-out_plot']=[x_mean*0.60,(x_out)*400-200,y_mean*0.6,(y_out+0.2)*200-100]
# 		elif x_out==0.0:srex_polygons[region]['line-out_plot']=[x_mean*0.60,(x_out+0.1)*400-200,y_mean*0.6,(y_out+0.1)*200-100]
# 		elif x_out==0.9:srex_polygons[region]['line-out_plot']=[x_mean*0.60,(x_out)*400-200,y_mean*0.6,(y_out+0.1)*200-100]
# 		elif y_out==0.0:srex_polygons[region]['line-out_plot']=[x_mean*0.60,(x_out+0.05)*400-200,y_mean*0.6,(y_out+0.2)*200-100]
# 		elif y_out==0.8:srex_polygons[region]['line-out_plot']=[x_mean*0.60,(x_out+0.05)*400-200,y_mean*0.6,(y_out)*200-100]



# for region in srex_polygons.keys():
# 	srex_polygons[region]['out_plot'][0]+=0.01
# 	srex_polygons[region]['out_plot'][1]+=0.005
# 	srex_polygons[region]['out_plot'][2]-=0.02
# 	srex_polygons[region]['out_plot'][3]-=0.06

# output = open('srex_dict.pkl', 'wb')
# pickle.dump(srex_polygons, output)
# output.close()

