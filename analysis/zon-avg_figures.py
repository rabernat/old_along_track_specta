from pylab import *
import h5py
import os
import mycolors
from scipy.ndimage.filters import gaussian_filter1d, gaussian_filter

# which data to use
prefix = { 'SAT': 'SAT_50degwide', 'POP': 'POP_50degwide' } 
# the different variables available
varnames = ['V','U','T','VT','VU']
cnames = ['k', 'om', 'c']
expts = ['SAT','POP']
# load data
data = dict()
dza = dict()
for ex in expts:
    data[ex] = dict()
    for v in varnames:
        data[ex][v] = dict(np.load('../data/%s_%s.npz' % (prefix[ex], v)))
    dza[ex] = dict(np.load('../data/%s_%s.npz' % (prefix[ex], 'zon-avg')))

# calculate moments
M1, M2, qual = dict(), dict(), dict()
for ex in expts:
    M1[ex], M2[ex], qual[ex] = dict(), dict(), dict()
    for v in varnames:
        M1[ex][v], M2[ex][v], qual[ex][v] = dict(), dict(), dict()
        for coord in cnames:
            thecoord = data[ex][v][coord]
            thedat = data[ex][v]['pow_' + coord].copy()
            thedat = ma.masked_array(thedat, zeros_like(thedat))
            # normalize
            normfac = sum(thedat, axis=1)
            # a check for profiles that have both positive and negative values
            qual[ex][v][coord] = (abs(normfac) / sum(abs(thedat), axis=1))
            if coord=='k':
                thedat.mask[:,:4] = True
            thedat /= normfac[:,newaxis]
            thedat = ma.masked_invalid(thedat)
            thedat.mask += (qual[ex][v][coord]<0.9)[:,newaxis]

            M1[ex][v][coord] = sum(thedat * thecoord, axis=1)
            M2[ex][v][coord] = (
                sum(thedat * (thecoord - M1[ex][v][coord][:,newaxis])**2, axis=1) )
            M2[ex][v][coord].mask += (M2[ex][v][coord] < 0.)

# load grid info from data
d = data['SAT']['T']
k, om, c = d['k'], d['om'], d['c']
dk, dom, dc = d['dk'], d['dom'], d['dc']
Nk, Nom, Nc = k.shape[1], len(om), c.shape[1]
lat = d['lat']
Nlat = len(lat)

# other datasets
andreas_data_dir = os.path.join(os.environ['D'], 'DATASTORE.RPA','projects','aviso_mixing','andreas')
cdat = np.load(os.path.join(andreas_data_dir, 'c.npz'))
clat = linspace(-80,80,160)
Udat = np.load(os.path.join(andreas_data_dir, 'Umean_ECCO_patch.npz'))
rdat = np.load(os.path.join(andreas_data_dir, 'r.npz'))

leg = [r'$\overline{|V|^2}$',r'$\overline{|\Theta|^2}$',r'$\overline{V\Theta}$']


rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 8
close('all')
figure(figsize=(6.5,7.75))
ax1=subplot(311)
plot( lat, 2*pi*M1['SAT']['V']['k']**-1 / 1e3, 'b-',
      lat, 2*pi*M1['SAT']['T']['k']**-1 / 1e3, 'g-',
      lat, 2*pi*M1['SAT']['VT']['k']**-1 / 1e3, 'r-',
      clat, rdat['r_dudley'] / 1e3, 'k-',
      clat, rdat['r_rossby'] / 1e3, 'k--',
      lat, 2*pi*M1['POP']['V']['k']**-1 / 1e3, 'b--',
      lat, 2*pi*M1['POP']['T']['k']**-1 / 1e3, 'g--',
      lat, 2*pi*M1['POP']['VT']['k']**-1 / 1e3, 'r--')
#xlabel('latitude')
ylabel(r'wavelength (km)')
xlim([-60,50])
ylim([0,2000])
legend(leg + [r'$L_{eddy}$', r'$L_{d}$'], loc='upper left')
title(r'$M_1^\kappa$')
grid()
ax2=subplot(312)
plot( lat, M2['SAT']['V']['k']**0.5 / (2*pi) * 1e3, 'b-',
      lat, M2['SAT']['T']['k']**0.5 / (2*pi) * 1e3, 'g-',
      lat, M2['SAT']['VT']['k']**0.5 / (2*pi) * 1e3, 'r-',
      lat, M2['POP']['V']['k']**0.5 / (2*pi) * 1e3, 'b--',
      lat, M2['POP']['T']['k']**0.5 / (2*pi) * 1e3, 'g--',
      lat, M2['POP']['VT']['k']**0.5 / (2*pi) * 1e3, 'r--', )
#xlabel('latitude')
ylabel(r'spectral width (cycles / km)')
xlim([-60,50])
#ylim([0,1000])
title(r'$\sqrt{M_2^\kappa}$')
grid()
# ratio between M1[k] and L_eddy
L_eddy = interp(lat, clat, rdat['r_dudley'])
L_ros = interp(lat, clat, rdat['r_rossby'])
subplot(313)
plot( lat, 2*pi*M1['SAT']['VT']['k']**-1 / L_eddy, 'k-',
      lat, 2*pi*M1['SAT']['VT']['k']**-1 / L_ros, 'c-',
      lat, 2*pi*M1['POP']['VT']['k']**-1 / L_eddy, 'k--',
      lat, 2*pi*M1['POP']['VT']['k']**-1 / L_ros, 'c--')
grid()
xlim([-60,50]); ylim([0,30])
xlabel('lat');
legend([r'$\gamma / L_{eddy}$', r'$\gamma / L_{d}$'],
        loc='upper center')
title('Ratio Between Length Scales')

tight_layout()
savefig('../figures/moments_k.pdf')

figure(figsize=(6.5,5.5))
ax1=subplot(211)
plot( lat, M1['SAT']['V']['c'], 'b-',
      lat, M1['SAT']['T']['c'], 'g-',
      lat, M1['SAT']['VT']['c'], 'r-',
      clat, -cdat['c_dudley'], 'k-',
      clat, cdat['c_doppler'], 'k--',
      lat, M1['POP']['V']['c'], 'b--',
      lat, M1['POP']['T']['c'], 'g--',
      lat, M1['POP']['VT']['c'], 'r--')
xlabel('latitude')
ylabel(r'phase speed (m s$^{-1}$)')
xlim([-60,50])
ylim([0.1,-0.5])
title(r'$M_1^c$')
legend(leg + [r'$c_{eddy}$', r'$c_{R}$'], loc='upper left')
grid()
ax2=subplot(212)
plot( lat, M2['SAT']['V']['c']**0.5, 'b-',
      lat, M2['SAT']['T']['c']**0.5, 'g-',
      lat, M2['SAT']['VT']['c']**0.5, 'r-',
      lat, M2['POP']['V']['c']**0.5, 'b--',
      lat, M2['POP']['T']['c']**0.5, 'g--',
      lat, M2['POP']['VT']['c']**0.5, 'r--' )
xlabel('latitude')
ylabel(r'phase speed (m s$^{-1}$)')
xlim([-60,50])
title(r'$\sqrt{M_2^c}$')
grid()
tight_layout()
savefig('../figures/moments_c.pdf')


# Holt & Talley MLD
mld_data_dir =  os.path.join(os.environ['D'], 'mixedlayer')
mld_hdf_file = h5py.File(os.path.join(mld_data_dir, 'climatology.nc'),'r')
# smooth data with gaussian filter
mld = gaussian_filter1d(
        ma.masked_invalid(mld_hdf_file.get('dt_mld_mean')[:,:30]).mean(axis=1), 2)
mld_lat = ma.masked_invalid(mld_hdf_file.get('latgrid')[:,0])
# a few other constants
rho0 = 1027.
cp = 4186.

mask = (sstmask | sshmask)
MHT = rho0*cp*s.L*interp(lat,mld_lat,mld) * ma.masked_array(
            za_data['VpTp'], mask) 
DY = (s.lat[1] - s.lat[0])*110e3
# 2nd order centered difference
MHT_smooth = MHT.filled(0.)
MHT_smooth[mask] = interp(lat[mask],s.lat[~mask],MHT[~mask])
MHT_smooth = gaussian_filter1d(MHT_smooth,1.5)
heating_rate = -hstack([0, (MHT_smooth[2:] - MHT_smooth[:-2]), 0]) / s.L / (2*DY)

Tbar = za_data['Tbar']
dTbar_dy = hstack([0, Tbar[2:] - Tbar[:-2] ,0]) / (2*DY)
K = ma.masked_array(-za_data['VpTp']/dTbar_dy, abs(dTbar_dy) < 1e-6)

 
#close('all')

# mht
figure(figsize=(6.5,5.5))
subplot(211)
plot(s.lat, MHT,'k',linewidth=2)
grid(); xlim([-60,50]); ylim(array([-1,1])*1.3e13)
xlabel('lat'); ylabel(r'MHT (W)')
title('Meridional Heat Transport')

# heating
subplot(212)
plot(s.lat, ma.masked_array(heating_rate, mask), 'k', linewidth=2)
grid(); xlim([-60,50]); ylim([-15,15])
xlabel('lat'); ylabel(r'Q (W m$^{-2}$)')
title('Heating Rate')
tight_layout()

savefig_dummy('../figures/%s/MHT_and_heating.pdf' % secname)

figure(figsize=(6.5,2.25))
plot(s.lat, K,'k',linewidth=2)
grid(); xlim([-60,50]); ylim(array([0,5000]))
xlabel('lat'); ylabel(r'K (m$^2$ s$^{-1}$)')
title('Meridional Eddy Diffusivity')
tight_layout()
savefig_dummy('../figures/%s/diffusivity.pdf' % secname)


# output Tbar for advection/diffusion calc
Tbar_fine = tile( interp(arange(-80,80,0.1)+0.5, s.lat, za_data['Tbar'])[:,newaxis],[1,500] )
#Tbar_fine.astype(dtype('>f4')).tofile('../data/PACE_SST.bin')
dTdy = (za_data['Tbar'][2:]-za_data['Tbar'][:-2]) / (2*DY)

alpha_c = - data['VT']['pow_c']/(data['V']['pow_c']**0.5 * data['T']['pow_c']**0.5)


