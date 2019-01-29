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
    

ccmp_path = '/home/orca/mcurcic/data/ccmp_on_merra2'

starttime = datetime(2016,12,1)
endtime = datetime(2017,3,1)

time = starttime
while True:

    time0, time1, w0, w1 = read_times(time) 

    if time == endtime:
        time1 = time

    with Dataset(ccmp_path+'/ccmp_on_merra2_'+time0.strftime('%Y-%m-%d_%H:%M:%S')+'.nc') as nc:
        u0 = nc.variables['u_ccmp'][:]
        v0 = nc.variables['v_ccmp'][:]

    with Dataset(ccmp_path+'/ccmp_on_merra2_'+time1.strftime('%Y-%m-%d_%H:%M:%S')+'.nc') as nc:
        u1 = nc.variables['u_ccmp'][:]
        v1 = nc.variables['v_ccmp'][:]

    jdm, idm = u0.shape

    u0[np.isnan(u0)] = 0
    v0[np.isnan(v0)] = 0
    u1[np.isnan(u1)] = 0
    v1[np.isnan(v1)] = 0

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
