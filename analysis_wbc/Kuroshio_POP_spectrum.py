from pylab import *
import os
from scipy.io import netcdf
import mpl_toolkits.basemap as bm
import pyresample

# the POP model data
popstr = 'hybrid_v5_rel04_BC5_ne120_t12_pop62'
data_dir = os.path.join('/Volumes/Bucket1/Data/', popstr, 'processed')
output_dir = os.path.join(data_dir, 'movie')
nc = netcdf.netcdf_file(os.path.join(data_dir, 'Kuroshio_Ext.nc'), version=2)

Nt,Ny,Nx = nc.variables['SST'].shape
DX = diff(nc.variables['x'][:])[0]


area_name = 'Kuroshio Orthographic Grid'
proj_id = 'ease_ortho'
area_id = proj_id
#proj4_args = '+proj=laea +lat_0=32 +lon_0=172 +a=6371228.0 +units=m'
proj4_args = '+proj=ortho +lat_0=32 +lon_0=172 +units=m'
x_size = 400
y_size = 400
area_extent = (-1600000.,-1600000.,1600000.,1600000.)
area_def = pyresample.utils.get_area_def(area_id, area_name, proj_id, proj4_args,
                                       x_size, y_size, area_extent)

k = fftshift(fftfreq(Nx, area_def.pixel_size_x))
l = fftshift(fftfreq(Ny, area_def.pixel_size_y))
dk = k[1] - k[0]
dl = l[1] - l[0]
R = (k[newaxis,:]**2 + l[:,newaxis]**2)**0.5 
K = k[Nx/2:] # coordinate for plotting, positive definite
Nk = len(K)
# for binning
Ridx = digitize(R.ravel(), K, right=True)

varnames = {'T': 'SST', 'U': 'U1_1', 'V': 'V1_1'}
power_spectrum = dict()
for vname in varnames.iterkeys():
    power_spectrum[vname] = zeros((Nt,Nk))
mld = zeros(Nt)

window = hanning(Nx)*hanning(Ny)[:,newaxis]

for n in range(Nt):
    print n
    for vname, ncvname in varnames.iteritems():
        mld[n] = nc.variables['HMXL_2'][n].mean()
        Q = nc.variables[ncvname][n].copy() / (dk*dl)
        # remove mean and window
        Q -= Q.mean()
        Q += window
        # fourier transform
        Qf = fftshift(fft2(Q))
        # isoptropic power spectrum
        power_spectrum[vname][n] = (
            bincount(Ridx, weights=real(Qf*conj(Qf)).ravel())[:-1] )
        
idx_winter = mod(T-2,12)<1
idx_summer = mod(T-8,12)<1

for vname in varnames.keys():
    power_spectrum[vname + '_winter'] = (
        power_spectrum[vname][idx_winter].mean(axis=0))
    power_spectrum[vname + '_summer'] = (
        power_spectrum[vname][idx_summer].mean(axis=0))

figure()
loglog(K[1:], power_spectrum['V_summer'][1:] + power_spectrum['U_summer'][1:])
loglog(K[1:], power_spectrum['V_winter'][1:] + power_spectrum['U_winter'][1:])
loglog(K[10:130], 1e17*K[10:130]**-3, 'k-')
loglog(K[10:130], 1e25*K[10:130]**-(5/3.), 'k--')

    