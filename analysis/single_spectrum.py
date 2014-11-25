import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import periodogram, welch
import netCDF4

#ddir = '/Users/rpa/RND/Data/hybrid_v5_rel04_BC5_ne120_t12_pop62/'
ddir = '/Volumes/Bucket1/Data/hybrid_v5_rel04_BC5_ne120_t12_pop62/'
fname = 'hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1.0046-01-01.nc'
nc = netCDF4.Dataset(ddir + fname)

n = 0
j = 1500
ir = np.r_[2500:3500]
T = nc.variables['SST'][n,j,ir]
dx = nc.variables['DXT'][j,ir]/100.
x = np.cumsum(dx)

SSTps_den = 0.
Vps_den = 0.
cnt = 0
for year in range(46,50):
    for month in range(1,13):
        fname = 'hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1.%04d-%02d-01.nc' % (year, month)
        print fname
        nc = netCDF4.Dataset(ddir + fname)
        for n in range(30):
            cnt += 1
            T = nc.variables['SST'][n,j,ir]
            V = nc.variables['V1_1'][n,j,ir]/100.
            k, tf = periodogram(T, dx.mean()**-1, window='hanning', detrend='linear')
            k, vf = periodogram(V, dx.mean()**-1, window='hanning', detrend='linear')
            SSTps_den += tf
            Vps_den += vf
SSTps_den /= cnt
Vps_den /= cnt

plt.figure(figsize=(10,4))
plt.subplot(121)
plt.loglog(k**-1, Vps_den)
#plt.loglog(k[-1], 1e3, 'ro')
plt.xlim([1e5, 1e4])
plt.ylim([1e-1, 1e4])
plt.xlabel('Wavelength (m)')
plt.title('Kinetic Energy Power Spectrum')

plt.subplot(122)
plt.loglog(k**-1, SSTps_den)
#plt.loglog(k[-1], 1e3, 'ro')
plt.xlim([1e5, 1e4])
plt.ylim([1e-1, 1e4])
plt.xlabel('Wavelength (m)')
plt.title('SST Power Spectrum')

plt.tight_layout()
plt.savefig('figures/POP_spectrum.pdf')