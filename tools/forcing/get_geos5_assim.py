#!/usr/bin/env python

import subprocess
from datetime import datetime,timedelta

starttime = datetime(2016,4, 9,0,30)
endtime = datetime(2016,4,13,23,30)
time = starttime

baseurl = 'ftp://gmao_ops@ftp.nccs.nasa.gov/fp/das'

while True:
    filename = 'GEOS.fp.asm.tavg1_2d_slv_Nx.'+time.strftime('%Y%m%d_%H%M')+'.V01.nc4'
    url = baseurl+'/'+time.strftime('Y%Y/M%m/D%d/')+filename
    print url
    subprocess.call('wget --password="" '+url,shell=True)
    time += timedelta(hours=1)
    if time > endtime:break
