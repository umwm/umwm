#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from datetime import datetime,timedelta

starttime = datetime(2016,2, 8,0)
endtime   = datetime(2016,2,11,0)
time = starttime

path = '/home/disk/uwin2/milan/laser/realtime/output/'+starttime.strftime('%Y%m%d%H')

while True:

    umwmout = path+'/umwmout_'+time.strftime('%Y-%m-%d_%H:%M:%S')+'.nc'
    print 'Processing '+umwmout

    nc = Dataset(umwmout,'r')
    lon = nc.variables['lon'][0,:,:]
    lat = nc.variables['lat'][0,:,:]
    uc = nc.variables['uc'][0,:,:]
    vc = nc.variables['vc'][0,:,:]
    wspd = nc.variables['wspd'][0,:,:]
    wdir = nc.variables['wdir'][0,:,:]
    rhoa = nc.variables['rhoa'][0,:,:]
    rhow = nc.variables['rhow'][0,:,:]
    uw = wspd*np.cos(wdir)
    vw = wspd*np.sin(wdir)
    nc.close()

    jdm,idm = uw.shape

    umwmin = 'umwmin_'+time.strftime('%Y-%m-%d_%H:%M:%S')+'.nc'

    nc = Dataset(umwmin,'w',format='NETCDF3_CLASSIC')
    nc.createDimension('x',size=idm)
    nc.createDimension('y',size=jdm)
    nc.createVariable('lon',datatype='f4',dimensions=('y','x'))[:,:] = lon[:,:]
    nc.createVariable('lat',datatype='f4',dimensions=('y','x'))[:,:] = lat[:,:]
    nc.createVariable('uw',datatype='f4',dimensions=('y','x'))[:,:] = uw[:,:]
    nc.createVariable('vw',datatype='f4',dimensions=('y','x'))[:,:] = vw[:,:]
    nc.createVariable('uc',datatype='f4',dimensions=('y','x'))[:,:] = uc[:,:]
    nc.createVariable('vc',datatype='f4',dimensions=('y','x'))[:,:] = vc[:,:]
    nc.createVariable('rhoa',datatype='f4',dimensions=('y','x'))[:,:] = rhoa[:,:]
    nc.createVariable('rhow',datatype='f4',dimensions=('y','x'))[:,:] = rhow[:,:]
    nc.close()

    time += timedelta(hours=1)
    if time > endtime:break
