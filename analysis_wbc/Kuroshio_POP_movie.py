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

dpi=80.
close('all')
rcParams['font.size'] = 14
figure(figsize=(1920/dpi,1080/dpi), dpi=dpi)
for n in range(Nt):
    SST = nc.variables['SST'][n]
    U = nc.variables['U1_1'][n]
    V = nc.variables['V1_1'][n]
    U_rms = (U**2 + V**2)**0.5
    dUdx,dUdy= gradient(U) / DX
    dVdx,dVdy= gradient(V) / DX
    Z = -dUdy + dVdx
    subplot(131)
    imshow(SST, cmap=get_cmap('RdYlBu_r'));
    clim([2,28]); xticks([]); yticks([]);
    title('SST')
    subplot(132)
    imshow(U_rms, cmap=get_cmap('CMRmap'))
    clim([0,1]); xticks([]); yticks([]);
    title('Current Speed')
    subplot(133)
    imshow(Z, cmap=get_cmap('bwr'))
    clim([-5e-6,5e-6]); xticks([]); yticks([]);
    title('Vorticity')
    tight_layout()
    savefig(os.path.join(output_dir,'Kuroshio_frame_%05d.png' % n),
            figsize=(1920/dpi,1080/dpi), dpi=dpi)
    
    