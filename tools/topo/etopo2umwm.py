#!/usr/bin/env python

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np

infile = 'umwm.grid'

etopofile = '/home/orca/mcurcic/data/etopo01/ETOPO1_Ice_c_gmt4.grd'
with Dataset(etopofile, 'r') as nc:
    elon = nc.variables['x'][:]
    elat = nc.variables['y'][:]
    ez = nc.variables['z'][:,:]

with Dataset(infile, 'r') as nc:
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]

dlon = np.diff(lon[0,:])[-1]
dlat = np.diff(lat[:,0])[-1]

delon = np.diff(elon)[-1]
delat = np.diff(elat)[-1]

dlon_fac = int(dlon / delon)
dlat_fac = int(dlat / delat)

print(dlon,dlat)
print(delon,delat)
print(dlon_fac,dlat_fac)

z = np.zeros((lon.shape))
jdm, idm = z.shape

lon1d = lon[0,:]
lat1d = lat[:,0]


jj = dlat_fac // 2
ii = dlon_fac // 2
print(ii,jj)

for j in range(jdm):
    j0 = np.argmin((lat1d[j]-elat)**2)
    print(j)
    for i in range(1,idm):
        i0 = np.argmin((lon1d[i]-elon)**2)
        z[j,i] = np.mean(ez[j0-jj:j0+jj,i0-ii:i0+ii])
    z[j,0] = np.mean(np.mean(ez[j0-jj:j0+jj,i0]))

z[340:,:] = 1000.

outfile = 'umwm.gridtopo'
print('Writing '+outfile)
with Dataset(outfile, 'w', format='NETCDF3_CLASSIC') as nc:
    nc.createDimension('x', size=idm)
    nc.createDimension('y', size=jdm)
    nc.createVariable('lon', datatype='f4', dimensions=('y', 'x'))[:] = lon[:,:]
    nc.createVariable('lat', datatype='f4', dimensions=('y', 'x'))[:] = lat[:,:]
    nc.createVariable('z', datatype='f4', dimensions=('y', 'x'))[:] = z[:,:]
