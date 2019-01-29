#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile',help='input MERRA-2 file')
args = parser.parse_args()

infile = args.infile

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
import os

starttime = datetime.strptime(os.path.basename(infile)[27:35], '%Y%m%d')

with Dataset(infile, 'r') as nc:
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    u = nc.variables['U10M'][:,:,:]
    v = nc.variables['V10M'][:,:,:]
    T = nc.variables['T2M'][:,:,:]
    q = nc.variables['QV2M'][:,:,:]
    p = nc.variables['SLP'][:,:,:]

lon, lat = np.meshgrid(lon, lat)

ndm, jdm, idm = u.shape

Rd = 286.9
Rw = 461.5

rhoa = p*(1 + q) / (T*(Rd + Rw*q))

for n in range(ndm):
    outfile = 'umwmin_'+(starttime+timedelta(hours=n)).strftime('%Y-%m-%d_%H:%M:%S')+'.nc'
    print('Writing '+outfile)
    with Dataset(outfile, 'w', format='NETCDF3_CLASSIC') as nc:
        nc.createDimension('x', size=idm)
        nc.createDimension('y', size=jdm)
        nc.createVariable('uw', datatype='f4', dimensions=('y', 'x'))[:] = u[n,:,:]
        nc.createVariable('vw', datatype='f4', dimensions=('y', 'x'))[:] = v[n,:,:]
        nc.createVariable('rhoa', datatype='f4', dimensions=('y', 'x'))[:] = rhoa[n,:,:]
