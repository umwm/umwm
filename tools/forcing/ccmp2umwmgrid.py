#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile',help='input CCMP file')
args = parser.parse_args()

infile = args.infile

from netCDF4 import Dataset
import numpy as np

with Dataset(infile, 'r') as nc:
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]

lon, lat = np.meshgrid(lon, lat)
jdm, idm = lon.shape

outfile = 'umwm.grid'
print('Writing '+outfile)
with Dataset(outfile, 'w', format='NETCDF3_CLASSIC') as nc:
    nc.createDimension('x', size=idm)
    nc.createDimension('y', size=jdm)
    nc.createVariable('lon', datatype='f4', dimensions=('y', 'x'))[:] = lon[:,:]
    nc.createVariable('lat', datatype='f4', dimensions=('y', 'x'))[:] = lat[:,:]
