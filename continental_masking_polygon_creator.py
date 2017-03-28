
import numpy as np
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import datetime

# IMPORT AND CONFIG 
import numpy as np
import netCDF4 as net
import dimarray as da 
import sys,glob,datetime,pickle,os,itertools
import pandas as pd
import matplotlib.pylab as plt 
from netCDF4 import Dataset,netcdftime,num2date

os.chdir('/Users/peterpfleiderer/Documents/')


m = Basemap()
m.readshapefile("masks/shapefiles/world/ne_50m_admin_0_countries", 'admin', drawbounds=False)

region_polygons={}

# Europe
name='Europe'
region_polygons.pop(name, None)
for shape, region in zip(m.admin, m.admin_info):
    region = {k.lower():v for k,v in region.items()}    
    if region['continent'] in ['Europe'] and region['admin']!='Russia':
        #print region['admin']
        if name in region_polygons.keys():
            region_polygons[name] = \
            region_polygons[name].symmetric_difference(Polygon(shape))
        else:
            region_polygons[name] = Polygon(shape)

region_polygons['Europe']=region_polygons['Europe'].intersection(Polygon([[0,35],[50,35],[50,70],[0,70]])).symmetric_difference(region_polygons['Europe'].intersection(Polygon([[-30,35],[0,35],[0,70],[-30,70]])))

# North_America
name='North_America'
region_polygons.pop(name, None)
for shape, region in zip(m.admin, m.admin_info):
    region = {k.lower():v for k,v in region.items()}    
    if region['admin'] in ['United States of America', 'Canada']:
        if name in region_polygons.keys():
            region_polygons[name] = \
            region_polygons[name].symmetric_difference(Polygon(shape))
        else:
            region_polygons[name] = Polygon(shape)

# Russia
name='Russia'
region_polygons.pop(name, None)
for shape, region in zip(m.admin, m.admin_info):
    region = {k.lower():v for k,v in region.items()}    
    if region['admin'] in ['Russia']:
        if name in region_polygons.keys():
            region_polygons[name] = \
            region_polygons[name].symmetric_difference(Polygon(shape))
        else:
            region_polygons[name] = Polygon(shape)

# Australia
name='Australia'
region_polygons.pop(name, None)
for shape, region in zip(m.admin, m.admin_info):
    region = {k.lower():v for k,v in region.items()}    
    if region['admin'] in ['Australia']:
        if name in region_polygons.keys():
            region_polygons[name] = \
            region_polygons[name].symmetric_difference(Polygon(shape))
        else:
            region_polygons[name] = Polygon(shape)

# Asia
name='Asia'
region_polygons.pop(name, None)
for shape, region in zip(m.admin, m.admin_info):
    region = {k.lower():v for k,v in region.items()}    
    if region['subregion'] in ['Eastern Asia', 'Southern Asia', 'Central Asia'] and region['admin']!='Russia':
        #print region['admin']
        #print region
        if name in region_polygons.keys():
            try:
                region_polygons[name] = \
                region_polygons[name].symmetric_difference(Polygon(shape))
            except:
                print 'problem with',region['admin']
        else:
            region_polygons[name] = Polygon(shape)

ax=Basemap()

ax.drawcoastlines()

try:
    for poly in region_polygons['Asia']:
        x, y = poly.exterior.xy
        ax.plot(x,y,'r-',c='g')
except:
    x, y = region_polygons['Asia'].exterior.xy
    ax.plot(x,y,'r-',c='g') 

plt.show()

# asdasdas


