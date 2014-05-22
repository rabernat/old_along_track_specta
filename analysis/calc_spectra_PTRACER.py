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
Nx = 500; Nxdata=None; secname = 'PSV_50degwide'

do_smoothing=True
sm_sig = 1.

s = sector_analyzer.Sector(Nx=Nx, lonmin=-179.95, dlon=0.1)
s.search_for_timeseries_data()
dsets = ['AVISO_passive_V-19930115_7day',
         'AVISO_passive_U-19930115_7day',
         'AVISO_passive_T-19930115_7day']

#Nc = 121
Nc = 1000
# we pick the phase speed grid now
cin = linspace(-1.,1.,Nc+1)
Nt = 486
Nk = s.Nk

# a way to bypass actually outputing anything
actually_save = True
def savefig_dummy(*args, **kwargs):
    if actually_save:
        savefig(*args, **kwargs)

# initialize results matrix
# fourier transform data
data = {'V':[],'U':[],'T':[],'VT':[],'VU':[]}
for v in data.keys():
    data[v] = {'pow_k': ma.masked_array(zeros((s.Ny, Nk)),True),
               'pow_om':  ma.masked_array(zeros((s.Ny, Nt)),True),
#               'pow_c':  ma.masked_array(zeros((s.Ny, Nc+2)),True),
#               'cpts':  ma.masked_array(zeros((s.Ny, Nc+2)),True) }
               'pow_c':  ma.masked_array(zeros((s.Ny, Nc)),True),
               'cpts':  ma.masked_array(zeros((s.Ny, Nc)),True) }
za_data = {'Vp2':zeros(s.Ny),'Tp2':zeros(s.Ny), 'Up2': zeros(s.Ny),
            'Tbar':zeros(s.Ny),'VpTp':zeros(s.Ny), 'VpUp':zeros(s.Ny)}

#c = zeros((s.Ny, Nc+2))
#dc = zeros((s.Ny, Nc+2))
c = zeros((s.Ny, Nc))
dc = zeros((s.Ny, Nc))

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
    
    savefig_dummy('../figures/%s/individual_spectra/SST_SST_wavefreq_spectra_%g.pdf' % (secname, int(round(SST.lat))) )
    
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
        data[v]['pow_c'][j], c[j], dc[j], data[v]['cpts'][j] = (
            SSH_V.sum_in_c_interp(field * mask, cin=cin, Nc=Nc) )
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
prefix = '../data/%s' % secname
for v in data.keys():
    np.savez('%s_%s.npz' % (prefix,v), lat=s.lat,
            c=c, dc=dc, k=s.k, dk=s.dk, om=SST.om, dom=SST.dom,
            pow_k=data[v]['pow_k'].filled(0.),
            pow_om=data[v]['pow_om'].filled(0.),
            pow_c=data[v]['pow_c'].filled(0.))
np.savez('%s_zon-avg.npz' % prefix, lat=s.lat, **za_data)
    
