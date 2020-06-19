#!/usr/bin/env python

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np

infile = 'umwm.grid'

topo_file = 'gebco_2019.nc'
with Dataset(topo_file, 'r') as nc:
    elon = nc.variables['lon'][:]
    elat = nc.variables['lat'][:]
    ez = nc.variables['elevation'][:,:]

with Dataset(infile, 'r') as nc:
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]

dlon = np.diff(lon[0,:])[-1]
dlat = np.diff(lat[:,0])[-1]

delon = np.diff(elon)[-1]
delat = np.diff(elat)[-1]

dlon_fac = int(dlon / delon)
dlat_fac = int(dlat / delat)

print(dlon, dlat)
print(delon, delat)
print(dlon_fac, dlat_fac)

z = np.zeros((lon.shape))
jdm, idm = z.shape

lon1d = lon[0,:]
lat1d = lat[:,0]


jj = dlat_fac // 2
ii = dlon_fac // 2
print(ii,jj)

for j in range(jdm):
    j0 = np.argmin((lat1d[j] - elat)**2)
    print('Processing row', j)
    for i in range(idm):
        i0 = np.argmin((lon1d[i]-elon)**2)
        z[j,i] = np.mean(ez[j0-jj:j0+jj,i0-ii:i0+ii])

outfile = 'umwm.gridtopo'
print('Writing '+outfile)
with Dataset(outfile, 'w', format='NETCDF3_CLASSIC') as nc:
    nc.createDimension('x', size=idm)
    nc.createDimension('y', size=jdm)
    nc.createVariable('lon', datatype='f4', dimensions=('y', 'x'))[:] = lon[:,:]
    nc.createVariable('lat', datatype='f4', dimensions=('y', 'x'))[:] = lat[:,:]
    nc.createVariable('z', datatype='f4', dimensions=('y', 'x'))[:] = z[:,:]
