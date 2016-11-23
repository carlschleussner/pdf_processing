#############
# SREX plot info
#############

pkl_file = open('srex_dict.pkl', 'rb')
srex_polygons = pickle.load(pkl_file)
pkl_file.close()    

srex_polygons['WNA']['out_plot']=[0.0,0.6,0.1,0.2]
srex_polygons['CAM']['out_plot']=[0.0,0.4,0.1,0.2]
srex_polygons['AMZ']['out_plot']=[0.0,0.2,0.1,0.2]

srex_polygons['WNA']['out_plot']=[0.0,0.8,0.1,0.2]
srex_polygons['ALA']['out_plot']=[0.125,0.8,0.1,0.2]
srex_polygons['CNA']['out_plot']=[0.25,0.8,0.1,0.2]
srex_polygons['CGI']['out_plot']=[0.375,0.8,0.1,0.2]
srex_polygons['NEU']['out_plot']=[0.5,0.8,0.1,0.2]
srex_polygons['CEU']['out_plot']=[0.625,0.8,0.1,0.2] 
srex_polygons['NAS']['out_plot']=[0.75,0.8,0.1,0.2]
srex_polygons['CAS']['out_plot']=[0.875,0.8,0.1,0.2]
	
srex_polygons['TIB']['out_plot']=[0.875,0.6,0.1,0.2]
srex_polygons['EAS']['out_plot']=[0.875,0.4,0.1,0.2]
srex_polygons['SEA']['out_plot']=[0.875,0.2,0.1,0.2]

srex_polygons['WSA']['out_plot']=[0.0,0.0,0.1,0.2]
srex_polygons['SSA']['out_plot']=[0.125,0.0,0.1,0.2]
srex_polygons['MED']['out_plot']=[0.25,0.0,0.1,0.2]
srex_polygons['SAH']['out_plot']=[0.375,0.0,0.1,0.2]
srex_polygons['WAS']['out_plot']=[0.5,0.0,0.1,0.2]
srex_polygons['SAS']['out_plot']=[0.625,0.0,0.1,0.2]
srex_polygons['SAU']['out_plot']=[0.75,0.0,0.1,0.2]
srex_polygons['NAU']['out_plot']=[0.875,0.0,0.1,0.2]



x_ratio=0.67
y_ratio=0.5

for region in srex_polygons.keys():
	if region not in ['EAF','SAF','WAF','global']:
		x,y=Polygon(srex_polygons[region]['points']).exterior.xy
		x_mean,y_mean=np.mean(x),np.mean(y)
		x_out,y_out=srex_polygons[region]['out_plot'][0],srex_polygons[region]['out_plot'][1]

		if x_out==0.0 and y_out==0.0: x_out+=0.125;y_out+=0.2
		elif x_out==0.875 and y_out==0.8: x_out+=0;y_out+=0
		elif x_out==0.875 and y_out==0.0: x_out+=0.0;y_out+=0.2
		elif x_out==0.0 and y_out==0.8: x_out+=0.125;y_out+=0

		elif x_out==0.0: x_out+=0.0625;y_out+=0.1
		elif x_out==0.875: x_out+=0.0;y_out+=0.1
		elif y_out==0.0: x_out+=0.0625;y_out+=0.2
		elif y_out==0.8: x_out+=0.0625;y_out+=0
		srex_polygons[region]['line-out_plot']=[x_mean*x_ratio,(x_out)*400-200,y_mean*y_ratio,(y_out)*200-100]


for region in srex_polygons.keys():
	srex_polygons[region]['out_plot'][0]+=0.01
	srex_polygons[region]['out_plot'][1]+=0.005
	srex_polygons[region]['out_plot'][2]+=0.005
	srex_polygons[region]['out_plot'][3]-=0.06

#del srex_polygons['global']

output = open('srex_dict.pkl', 'wb')
pickle.dump(srex_polygons, output)
output.close()



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

