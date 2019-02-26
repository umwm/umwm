#!/usr/bin/env python

import argparse

from netCDF4 import Dataset
import numpy as np


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Creates 'umwm.grid' file.")
    parser.add_argument('file', help='input file')
    args = parser.parse_args()
    
    infile = args.infile
    
    with Dataset(infile, 'r') as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
    
    lon, lat = np.meshgrid(lon, lat)
    jdm, idm = lon.shape
    
    outfile = 'umwm.grid'
    print('Writing '+outfile)
    with Dataset(outfile, 'w', format='NETCDF3_CLASSIC') as nc:
        nc.createDimension('x', size=idm)
        nc.createDimension('y', size=jdm)
        nc.createVariable('lon', datatype='f4', dimensions=('y', 'x'))[:] = lon[:,:]
        nc.createVariable('lat', datatype='f4', dimensions=('y', 'x'))[:] = lat[:,:]

