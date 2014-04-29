from pylab import *
import mycolors
import sector_analyzer
import os
import h5py
from scipy.ndimage.filters import gaussian_filter1d, gaussian_filter

# the narrow sector
#Nx = 120; secname = '30degwide'
# should be the same
#Nx = 120; Nxdata=200; secname='30degwide'
# the wide sector
Nx = 200; Nxdata=None; secname = '50degwide'

do_smoothing=True
sm_sig = 1.

s = sector_analyzer.Sector(Nx=Nx, Nxdata=Nxdata)
s.search_for_timeseries_data()
dsets = ['AVISO_dt_ref_global_merged_msla_v-20020605_7day',
         'AVISO_dt_ref_global_merged_msla_u-20020605_7day',
         'NCDC_AVHRR_AMSR_OI-20020605_7day']

Nc = 121
Nt = 486
Nk = s.Nk

# a way to bypass actually outputing anything
actually_save = True
def savefig_dummy(*args, **kwargs):
    if actually_save:
        savefig(*args, **kwargs)


# other datasets
andreas_data_dir = os.path.join(os.environ['D'], 'DATASTORE.RPA','projects','aviso_mixing','andreas')
cdat = np.load(os.path.join(andreas_data_dir, 'c.npz'))
clat = linspace(-80,80,160)
Udat = np.load(os.path.join(andreas_data_dir, 'Umean_ECCO_patch.npz'))
rdat = np.load(os.path.join(andreas_data_dir, 'r.npz'))

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

# initialize results matrix
# fourier transform data
data = {'V':[],'U':[],'T':[],'VT':[],'VU':[]}
for v in data.keys():
    data[v] = {'pow_k': ma.masked_array(zeros((s.Ny, Nk)),True),
               'pow_om':  ma.masked_array(zeros((s.Ny, Nt)),True),
               'pow_c':  ma.masked_array(zeros((s.Ny, Nc+2)),True),
               'cpts':  ma.masked_array(zeros((s.Ny, Nc+2)),True) }
za_data = {'Vp2':zeros(s.Ny),'Tp2':zeros(s.Ny), 'Up2': zeros(s.Ny),
            'Tbar':zeros(s.Ny),'VpTp':zeros(s.Ny), 'VpUp':zeros(s.Ny)}

c = zeros((s.Ny, Nc+2))
dc = zeros((s.Ny, Nc+2))

mask = ones(s.Nk)
mask[0] = 0

# figure stuff
day = 24*60*60.
#days = array([-15,-30,-60,-120,120,60,30,15])
days = array([-15,-30,-60,inf,60,30,15])
omtick = (day * days.astype('f4') / 2 / pi)**-1
lens =  array([1000,350,200,100,80,50])
ktick = (1000. * lens.astype('f4') / 2 / pi)**-1

rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 'small'
rcParams['axes.formatter.limits'] = [-2, 2]
def spectral_plot(SST,SSH_V):
    fig=figure(figsize=(6.5,7.5))

    subplot(311)
    pcolormesh(SSH_V.om, SSH_V.k, log10(real(SSH_V.ft_data*conj(SSH_V.ft_data)).T), rasterized=True)
    clim([-9,-5])
    xticks(omtick,days)
    yticks(ktick,lens)
    ylim([0,6.5e-5]); xlim(array([-1,1])*5e-6)
    xlabel(r'$2 \pi / \omega$ (days)'); ylabel(r'$2 \pi / \kappa)$ (km)')
    title(r'log$_{10}$ $|V|^2(\kappa,\omega)$ (m$^2$ s$^{-2}$)')
    colorbar()
    
    subplot(312)
    pcolormesh(SST.om, SST.k, log10(real(SST.ft_data*conj(SST.ft_data)).T), rasterized=True)
    xticks(omtick,days)
    yticks(ktick,lens)
    ylim([0,6.5e-5]); xlim(array([-1,1])*5e-6)
    clim([-7,-3])
    xlabel(r'$2 \pi / \omega$ (days)'); ylabel(r'$2 \pi / \kappa)$ (km)')
    title(r'log$_{10}$ $|\Theta|^2(\kappa,\omega)$ (K$^2$)')
    colorbar()
    
    subplot(313)
    pcolormesh(SST.om, SST.k, real(SSH_V.ft_data*conj(SST.ft_data)).T, 
        cmap=get_cmap('bwr'), rasterized=True)
    xticks(omtick,days)
    yticks(ktick,lens)
    ylim([0,6.5e-5]); xlim(array([-1,1])*5e-6)
    clim(array([-1,1])*1e-5)
    xlabel(r'$2 \pi / \omega$ (days)'); ylabel(r'$2 \pi / \kappa)$ (km)')
    title(r'Re$(V \Theta^\ast)$ (K m s$^{-1}$)')
    colorbar()
    
    fig.tight_layout()
    
    savefig_dummy('../figures/SAT_%s/individual_spectra/SST_SST_wavefreq_spectra_%g.pdf' % (secname, int(round(SST.lat))) )
    
plot_js = arange(79,s.Ny,40)
#plot_js = array([])

sstmask = zeros(s.Ny,bool)
sshmask = zeros(s.Ny,bool)
for j in arange(s.Ny):
    print(j)
    
    SSH_V = s.timeseries_sets[dsets[0]].load(j, remove_zonal_mean=True, Nt=Nt, ft_normfac=sqrt(2)/(Nt*s.Nx))
    SSH_U = s.timeseries_sets[dsets[1]].load(j, remove_zonal_mean=True, Nt=Nt, ft_normfac=sqrt(2)/(Nt*s.Nx))
    SST = s.timeseries_sets[dsets[2]].load(j, Nt=Nt, 
        remove_temporal_mean=True, remove_zonal_mean=True, ft_normfac=sqrt(2)/(Nt*s.Nx))
    
    if any(plot_js==j):
        spectral_plot(SST,SSH_V)
    
    myfields = []
    if any(SST.ts_data.sum(axis=0)==0.):
        sstmask[j] = True
    else:
        myfields.append(('T',  real(SST.ft_data*SST.ft_data.conj()) ))
    if any(SSH_V.ts_data.sum(axis=0)>1e10):
        sshmask[j] = True
    else:
        myfields.append(('V',  real(SSH_V.ft_data*SSH_V.ft_data.conj()) ))
        myfields.append(('U',  real(SSH_U.ft_data*SSH_U.ft_data.conj()) ))
        myfields.append(('VU',  real(SSH_U.ft_data*SSH_V.ft_data.conj()) ))
    if (not sstmask[j]) and (not sshmask[j]):
        myfields.append(('VT', 2*real(SSH_V.ft_data*SST.ft_data.conj())))
    VTf = 2 * real( SSH_V.ft_data * SST.ft_data.conj() )
    for (v, field) in myfields:
        if do_smoothing:
            field = gaussian_filter(field, sm_sig)
        data[v]['pow_k'][j] = SSH_V.sum_over_om(field * mask)
        data[v]['pow_om'][j] = SSH_V.sum_over_k(field * mask)
        #data[v]['pow_c'][j], c[j], dc[j], data[v]['cpts'][j] = SSH_V.sum_in_c(field * mask, Nc)
        data[v]['pow_c'][j], c[j], dc[j], data[v]['cpts'][j] = SSH_V.sum_in_c_interp(field * mask, Nc=Nc)
    # zero wavenumber and frequency are already removed
    #Vp = SSH_V.ts_data - SSH_V.ts_data.mean(axis=1)[:,newaxis]      
    #Up = SSH_U.ts_data - SSH_U.ts_data.mean(axis=1)[:,newaxis]      
    #Tp = SST.ts_data - SST.ts_data.mean(axis=1)[:,newaxis]
    # zonal and time mean
    Tbar = ma.masked_equal(SST.ts_data,0.).mean()
    za_data['Tbar'][j] = Tbar
    za_data['Vp2'][j] = mean(SSH_V.ts_data_filtered**2)
    za_data['Up2'][j] = mean(SSH_U.ts_data_filtered**2)
    za_data['Tp2'][j] = mean(SST.ts_data_filtered**2)
    za_data['VpTp'][j] = mean(SSH_V.ts_data_filtered*SST.ts_data_filtered)
    za_data['VpUp'][j] = mean(SSH_V.ts_data_filtered*SSH_U.ts_data_filtered)

# save data
prefix = '../data/SAT_%s' % secname
for v in data.keys():
    np.savez('%s_%s.npz' % (prefix,v), lat=s.lat,
            c=c, dc=dc, k=s.k, dk=s.dk, om=SST.om, dom=SST.dom,
            pow_k=data[v]['pow_k'].filled(0.),
            pow_om=data[v]['pow_om'].filled(0.),
            pow_c=data[v]['pow_c'].filled(0.))
np.savez('%s_zon-avg.npz' % prefix, lat=s.lat, **za_data)
    

mask = (sstmask | sshmask)
MHT = rho0*cp*s.L*interp(s.lat,mld_lat,mld) * ma.masked_array(
            za_data['VpTp'], mask) 
DY = (s.lat[1] - s.lat[0])*110e3
# 2nd order centered difference
MHT_smooth = MHT.filled(0.)
MHT_smooth[mask] = interp(s.lat[mask],s.lat[~mask],MHT[~mask])
MHT_smooth = gaussian_filter1d(MHT_smooth,1.5)
heating_rate = -hstack([0, (MHT_smooth[2:] - MHT_smooth[:-2]), 0]) / s.L / (2*DY)

Tbar = za_data['Tbar']
dTbar_dy = hstack([0, Tbar[2:] - Tbar[:-2] ,0]) / (2*DY)
K = ma.masked_array(-za_data['VpTp']/dTbar_dy, abs(dTbar_dy) < 1e-6)

 
close('all')

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







    
