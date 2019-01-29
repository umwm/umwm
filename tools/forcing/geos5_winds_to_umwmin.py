#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
from datetime import datetime,timedelta

starttime = datetime(2016,4, 9,0,30)
endtime = datetime(2016,4,13,23,30)
time = starttime

gridfile = 'GEOS.fp.asm.tavg1_2d_slv_Nx.20160409_0030.V01.nc4'

nc = Dataset(gridfile)
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
nc.close()

lon,lat = np.meshgrid(lon,lat)

jdm,idm = lon.shape

while True:

    file = 'GEOS.fp.asm.tavg1_2d_slv_Nx.'+time.strftime('%Y%m%d_%H%M')+'.V01.nc4'
    print file
  
    with Dataset(file,'r') as nc:
        u = nc.variables['U2M'][0,:,:]
        v = nc.variables['V2M'][0,:,:]

    print time,np.mean(u),np.mean(v)

    nc = Dataset('umwmin_'+time.strftime('%Y-%m-%d_%H:00:%S')+'.nc',\
                 'w',format='NETCDF3_CLASSIC')
    nc.createDimension('x',size=idm)
    nc.createDimension('y',size=jdm)
    nc.createVariable('uw',datatype='f4',dimensions=('y','x'))[:] = u[:]
    nc.createVariable('vw',datatype='f4',dimensions=('y','x'))[:] = v[:]
    nc.close()

    time += timedelta(hours=1)
    if time > endtime:break


