#!/usr/bin/env python

from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np

infile = '/gpfsm/dnb32/araman2/sandbox/experiments/offline_umwm_seaice/umwm_seaice_new/umwm_work/src/input/umwm.grid'

transparencyfile = '/gpfsm/dnb32/araman2/sandbox/experiments/offline_umwm_seaice/umwm_seaice_new/umwm_work/src/input/umwm_transparency4.nc'

nc   = Dataset(transparencyfile, 'r') 
glon = nc.variables['longitude'][:]
glat = nc.variables['latitude'][:]
gz   = nc.variables['depths'][:,:,:]
gsx  = nc.variables['sx1'][:,:,:]
gsy  = nc.variables['sy1'][:,:,:]
gz = np.squeeze(gz)
gsx = np.squeeze(gsx)
gsy = np.squeeze(gsy)
nc.close()

nc = Dataset(infile, 'r')
lon = nc.variables['lon'][:,:]
lat = nc.variables['lat'][:,:]
nc.close()

dlon = np.diff(lon[0,:])[-1]
dlat = np.diff(lat[:,0])[-1]

delon = np.diff(glon)[-1]
delat = np.diff(glat)[-1]

dlon_fac = int(dlon / delon)
dlat_fac = int(dlat / delat)

print(dlon,dlat)
print(delon,delat)
print(dlon_fac,dlat_fac)

z  = np.zeros((lon.shape))
sx = np.zeros((lon.shape))
sy = np.zeros((lon.shape))

jdm, idm = z.shape

lon1d = lon[0,:]
lat1d = lat[:,0]


jj = dlat_fac // 2
ii = dlon_fac // 2
print(ii,jj)

for j in range(jdm):
    j0 = np.argmin((lat1d[j]-glat)**2)
    print(j)
    for i in range(1,idm):
        i0 = np.argmin((lon1d[i]-glon)**2)
        z[j,i]  = np.mean(gz[j0,i0])
        sx[j,i] = np.mean(gsx[j0,i0])
        sy[j,i] = np.mean(gsy[j0,i0])
    z[j,0]  = np.mean(np.mean(gz[j0,i0]))
    sx[j,0] = np.mean(np.mean(gsx[j0,i0]))
    sy[j,0] = np.mean(np.mean(gsy[j0,i0]))

z[340:,:] = 1000.

outfile = 'umwm.gridtopo'

print('Writing '+outfile)

nc = Dataset(outfile, 'w', format='NETCDF4_CLASSIC') 
nc.createDimension('x', size=idm)
nc.createDimension('y', size=jdm)
nc.createVariable('lon', datatype='f4', dimensions=('y', 'x'))[:] = lon[:,:]
nc.createVariable('lat', datatype='f4', dimensions=('y', 'x'))[:] = lat[:,:]
nc.createVariable('z', datatype='f4', dimensions=('y', 'x'))[:] = z[:,:]
nc.close()



outfile = 'umwm.gridtransparency'

print('Writing '+outfile)

nc = Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
nc.createDimension('x', size=idm)
nc.createDimension('y', size=jdm)
nc.createVariable('lon', datatype='f4', dimensions=('y', 'x'))[:] = lon[:,:]
nc.createVariable('lat', datatype='f4', dimensions=('y', 'x'))[:] = lat[:,:]
nc.createVariable('sx', datatype='f4', dimensions=('y', 'x'))[:] = sx[:,:]
nc.createVariable('sy', datatype='f4', dimensions=('y', 'x'))[:] = sy[:,:]
nc.close()
