from pylab import *
import os
from scipy.io import netcdf
from scipy.interpolate import interp2d
import pyresample

N = 60
dTdt = zeros(N) # rms change in sst
t = zeros(N)

popstr = 'hybrid_v5_rel04_BC5_ne120_t12_pop62'
data_dir = os.path.join('/Volumes/Bucket1/Data/', popstr)
# a sample netcdf file for loading grid information
ncfname = [popstr + '.pop.h.nday1.0049-09-01.nc',
           popstr + '.pop.h.nday1.0049-10-01.nc']
n = 0
nplot = 31

figure()
for fname in ncfname:
    nc = netcdf.netcdf_file(os.path.join(data_dir, fname))
    for nt in range(30):
        t[n] = nc.variables['time'][nt]
        SST = ma.masked_greater(nc.variables['SST'][nt],1e30)
        SST.mask += (SST==-1.0)
        if n > 0:
            dSST = SST-SST_prev
            dt = (t[n] - t[n-1])
            dTdt[n] = mean(dSST**2)**0.5 / dt
            print 'dTdt = %g' % dTdt[n]
        SST_prev = SST.copy() 
        if (n>=nplot-2) and (n<=nplot):
            subplot(1,3,n-nplot, axisbg='0.5')
            imshow(dSST/dt, origin='bottom', cmap=get_cmap('bwr'))
            clim([-1,1])
            colorbar(orientation='horizontal', ticks=[-1,0,1],shrink=0.5)
            xticks([]); yticks([])
            title('dT/dt (deg. per day) : time = %g' % t[n])
            
        
        n += 1
        
