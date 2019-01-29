#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile',help='input CCMP file')
args = parser.parse_args()

infile = args.infile

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
import os

starttime = datetime.strptime(os.path.basename(infile)[19:27], '%Y%m%d')

with Dataset(infile, 'r') as nc:
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]
    u = nc.variables['uwnd'][:,:,:]
    v = nc.variables['vwnd'][:,:,:]

lon, lat = np.meshgrid(lon, lat)

ndm, jdm, idm = u.shape

for n in range(ndm):
    outfile = 'umwmin_'+(starttime+timedelta(hours=6*n)).strftime('%Y-%m-%d_%H:%M:%S')+'.nc'
    print('Writing '+outfile)
    with Dataset(outfile, 'w', format='NETCDF3_CLASSIC') as nc:
        nc.createDimension('x', size=idm)
        nc.createDimension('y', size=jdm)
        nc.createVariable('uw', datatype='f4', dimensions=('y', 'x'))[:] = u[n,:,:]
        nc.createVariable('vw', datatype='f4', dimensions=('y', 'x'))[:] = v[n,:,:]
