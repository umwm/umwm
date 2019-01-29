#!/usr/bin/env python

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
import os

def read_times(time):
    """Given input time with hour, returns the closest start and end times 
    of a 6-hour interval and interpolation weights."""
    time0 = datetime(time.year, time.month, time.day, (time.hour // 6)*6)
    time1 = time0 + timedelta(hours=6)
    w1 = (time-time0).total_seconds()/(6*3600)
    w0 = 1-w1
    return time0, time1, w0, w1
    

ccmp_path = '/home/orca/mcurcic/data/ccmp'

starttime = datetime(2016,12,1)
endtime = datetime(2017,3,1)

time = starttime
while True:

    time0, time1, w0, w1 = read_times(time) 

    with Dataset(ccmp_path+'/CCMP_Wind_Analysis_'+time0.strftime('%Y%m%d')+'_V02.0_L3.0_RSS.nc') as nc:
        u0 = nc.variables['uwnd'][time0.hour // 6,:,:]
        v0 = nc.variables['vwnd'][time0.hour // 6,:,:]

    with Dataset(ccmp_path+'/CCMP_Wind_Analysis_'+time1.strftime('%Y%m%d')+'_V02.0_L3.0_RSS.nc') as nc:
        u1 = nc.variables['uwnd'][time1.hour // 6,:,:]
        v1 = nc.variables['vwnd'][time1.hour // 6,:,:]

    jdm, idm = u0.shape

    #print(time, time0, time1, w0, w1)
    print(time)

    outfile = 'umwmin_'+time.strftime('%Y-%m-%d_%H:%M:%S')+'.nc'
    with Dataset(outfile, 'w', format='NETCDF3_CLASSIC') as nc:
        nc.createDimension('x', size=idm)
        nc.createDimension('y', size=jdm)
        nc.createVariable('uw', datatype='f4', dimensions=('y', 'x'))[:] = u0*w0 + u1*w1
        nc.createVariable('vw', datatype='f4', dimensions=('y', 'x'))[:] = v0*w0 + v1*w1

    time += timedelta(hours=1)
    if time > endtime:break
