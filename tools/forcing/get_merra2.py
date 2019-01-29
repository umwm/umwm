#!/usr/bin/env python

from datetime import datetime,timedelta
import os
import subprocess

username = os.environ['DISC_USERNAME']
password = os.environ['DISC_PASSWORD']

starttime = datetime(2016,12,1)
endtime = datetime(2017,3,1)

baseurl = 'http://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I1NXASM.5.12.4'

vars = ['U10M','V10M','SLP','T2M','QV2M','time','lon','lat']

time = starttime

while True:

    filename = 'MERRA2_400.inst1_2d_asm_Nx.'+time.strftime('%Y%m%d')+'.nc'
    url = baseurl+'/'+time.strftime('%Y/%m')+'/MERRA2_400.inst1_2d_asm_Nx.'+time.strftime('%Y%m%d')+'.nc4.nc?'+','.join(vars)

    subprocess.call('wget '+url+' -O '+filename+' --user '+username+' --password '+password,shell=True)

    time += timedelta(days=1)
    if time > endtime:break
