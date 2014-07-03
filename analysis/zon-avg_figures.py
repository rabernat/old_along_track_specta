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
        try:
            data[ex][v] = dict(np.load('../data/%s_%s.npz' % (prefix[ex], v)))
        except KeyError:
            pass
    dza[ex] = dict(np.load('../data/%s_%s.npz' % (prefix[ex], 'zon-avg')))

# calculate diffusivities
K = dict()
sig=4
for ex in expts:
    K[ex] = dict()
    for tr in ['T','S']:
        try:
            vq = dza[ex]['Vp%sp' % tr]
            # mask extreme values, don't know how they got there
            vq = ma.masked_array(vq, abs(vq) > 1e5)
            vq = ma.masked_array( gaussian_filter1d(vq.filled(0.), sig), vq.mask)[1:-1]
            qbar = gaussian_filter1d(dza[ex]['%sbar' % tr], sig)
            lat = dza[ex]['lat']
            # centered difference
            dqdy = (qbar[2:] - qbar[:-2]) / (110e3 * (lat[2:] - lat[:-2]))
            # mask for diffusivity
            dmask = dqdy**2 < (dqdy**2).mean()/100.
            K[ex][tr] = {'flux': vq,
                         'grad': dqdy,
                         'diff': -ma.masked_array(vq / dqdy, dmask)}
        except KeyError:
            pass

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
K_ka = np.load(os.path.join(andreas_data_dir, 'K.npz'))
EKEdat = np.load(os.path.join(andreas_data_dir,'aviso_EKE.npz'))
u_rms = gaussian_filter1d(((EKEdat['U2mean'] + EKEdat['V2mean'])[540:690]**0.5).mean(axis=0),2)

# Earth stuff
Om = 7.292e-5
L = 6.371e6
Beta = 2*Om*cos(pi*lat/180.)/L

Kdef = rdat['r_rossby']**-1
# obs scale
#Kobs = ((2*rdat['r_dudley'] )/ (2*pi))**-1
#Kobs = 0.5*(rdat['r_dudley'])**-1
Kobs = (rdat['r_dudley']*sqrt(2))**-1
# for Rhines scale

#L_rhines = (interp(lat, EKEdat['lat'], u_rms)/Beta)**0.5
#Krhines = (L_rhines/ (2*pi))**-1
Krhines = sqrt( Beta / (2*interp(lat, EKEdat['lat'], u_rms)) )

L_eddy = interp(lat, clat, 2*pi/Kobs)
L_ros = interp(lat, clat, 2*pi/Kdef)
L_rhines =  2*pi/Krhines

leg = [r'$\overline{|V|^2}$',r'$\overline{|\Theta|^2}$',r'$\overline{V\Theta}$']


rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 8
close('all')
figure(figsize=(6.5,7.75))
ax1=subplot(311)
plot( lat, 2*pi*M1['SAT']['V']['k']**-1 / 1e3, 'b-',
      lat, 2*pi*M1['SAT']['T']['k']**-1 / 1e3, 'g-',
      lat, 2*pi*M1['SAT']['VT']['k']**-1 / 1e3, 'r-',
      lat, L_eddy / 1e3, 'k-',
      lat, L_ros / 1e3, 'k--',
      lat, 2*pi*M1['POP']['V']['k']**-1 / 1e3, 'b--',
      lat, 2*pi*M1['POP']['T']['k']**-1 / 1e3, 'g--',
      lat, 2*pi*M1['POP']['VT']['k']**-1 / 1e3, 'r--')
#xlabel('latitude')
ylabel(r'wavelength (km)')
xlim([-60,50])
ylim([0,2000])
legend(leg + [r'$L_{eddy}$', r'$L_{d}$'], loc='upper left')
title(r'$2 \pi / M_1^k$')
grid()
ax2=subplot(312)
plot( lat, M2['SAT']['V']['k']**0.5 / (2*pi) * 1e6, 'b-',
      lat, M2['SAT']['T']['k']**0.5 / (2*pi) * 1e6, 'g-',
      lat, M2['SAT']['VT']['k']**0.5 / (2*pi) * 1e6, 'r-',
      lat, M2['POP']['V']['k']**0.5 / (2*pi) * 1e6, 'b--',
      lat, M2['POP']['T']['k']**0.5 / (2*pi) * 1e6, 'g--',
      lat, M2['POP']['VT']['k']**0.5 / (2*pi) * 1e6, 'r--', )
#xlabel('latitude')
ylabel(r'spectral width (cycles / 1000 km)')
xlim([-60,50])
#ylim([0,1000])
title(r'$\sqrt{M_2^k} / 2 \pi$')
grid()
# ratio between M1[k] and L_eddy

subplot(313)
plot( lat, 2*pi*M1['SAT']['VT']['k']**-1 / L_eddy, 'k-',
      lat, 2*pi*M1['SAT']['VT']['k']**-1 / L_ros, 'c-',
      lat, 2*pi*M1['SAT']['VT']['k']**-1 / L_rhines, 'm-',
      lat, 2*pi*M1['POP']['VT']['k']**-1 / L_eddy, 'k--',
      lat, 2*pi*M1['POP']['VT']['k']**-1 / L_ros, 'c--',
      lat, 2*pi*M1['POP']['VT']['k']**-1 / L_rhines, 'm--',)
grid()
xlim([-60,50]); ylim([0,4])
xlabel('lat');
legend([r'$K_{eddy}/M_1^k(\overline{V\Theta})$',
        r'$K_{d}/M_1^k(\overline{V\Theta})$',
        r'$K_{\beta}/M_1^k(\overline{V\Theta})$'],
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
ylabel(r'spectral width (m s$^{-1}$)')
xlim([-60,50])
title(r'$\sqrt{M_2^c}$')
grid()
tight_layout()
savefig('../figures/moments_c.pdf')


# Holt & Talley MLD
mld_data_dir =  os.path.join(os.environ['D'], 'mixedlayer')
mld_hdf_file = h5py.File(os.path.join(mld_data_dir, 'climatology.nc'),'r')
# smooth data with gaussian filter
mld_ht = gaussian_filter1d(
        ma.masked_invalid(mld_hdf_file.get('dt_mld_mean')[:,:30]).mean(axis=1), 2)
mld_lat = ma.masked_invalid(mld_hdf_file.get('latgrid')[:,0])

# a few other constants
rho0 = 1027.
cp = 4186.
mld = 50.
#L = 50*110e3*cos(pi*lat[1:-1]/180.)
L = 110e3*cos(pi*lat[1:-1]/180.) # one degree
Qfac = rho0*cp*mld*L
 
#close('all')

# mht
figure(figsize=(6.5,4.5))

subplot(221)
plot(lat[1:-1], Qfac*K['SAT']['T']['flux']/1e12,'k',
     lat[1:-1], Qfac*K['POP']['T']['flux']/1e12,'k--')
grid(); xlim([-60,50]); ylim(array([-1,1]))
legend(['SAT','POP'], loc='upper left')
ylabel(r'(TW / deg. long.)')
title('Meridional Heat Transport')

# gradient
subplot(222)
plot(lat[1:-1], K['SAT']['T']['grad']*1e3,'k',
     lat[1:-1], K['POP']['T']['grad']*1e3,'k--')
grid(); xlim([-60,50]);
ylabel(r'$\partial \overline{T} / \partial \varphi$ ($^\circ$ / km)')
title('Mean Meridional SST Gradient')

# EKE
EKE_SAT = 0.5*( dza['SAT']['Vp2'] + dza['SAT']['Up2'])
EKE_POP = 0.5*( dza['POP']['Vp2'] + dza['POP']['Up2'])
subplot(223)
plot(lat, EKE_SAT,'k',
     lat, EKE_POP,'k--')
grid(); xlim([-60,50]);
xlabel('lat'); ylabel(r'(m$^2$ s$^{-2}$)')
title('Eddy Kinetic Energy')

# diffusivity
subplot(224)
plot(lat[1:-1], K['SAT']['T']['diff'],'k',
     lat[1:-1], K['POP']['T']['diff'],'k--')
grid(); xlim([-60,50]); ylim([-0,7000])
xlabel('lat'); ylabel(r'(m$^2$ s$^{-1}$)')
title('Eddy Diffusivity')

tight_layout()
savefig('../figures/MHT_gradient_EKE_diffusivity.pdf')

show()




