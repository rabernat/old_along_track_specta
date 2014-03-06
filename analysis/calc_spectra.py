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
    
    savefig('../figures/%s/individual_spectra/SST_SST_wavefreq_spectra_%g.pdf' % (secname, int(round(SST.lat))) )
    
#plot_js = arange(79,s.Ny,40)
plot_js = array([])

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
        data[v]['pow_c'][j], c[j], dc[j], data[v]['cpts'][j] = SSH_V.sum_in_c(field * mask, Nc)
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

savefig('../figures/%s/MHT_and_heating.pdf' % secname)

figure(figsize=(6.5,2.25))
plot(s.lat, K,'k',linewidth=2)
grid(); xlim([-60,50]); ylim(array([0,5000]))
xlabel('lat'); ylabel(r'K (m$^2$ s$^{-1}$)')
title('Meridional Eddy Diffusivity')
tight_layout()
savefig('../figures/%s/diffusivity.pdf' % secname)


# output Tbar for advection/diffusion calc
Tbar_fine = tile( interp(arange(-80,80,0.1)+0.5, s.lat, za_data['Tbar'])[:,newaxis],[1,500] )
#Tbar_fine.astype(dtype('>f4')).tofile('../data/PACE_SST.bin')
dTdy = (za_data['Tbar'][2:]-za_data['Tbar'][:-2]) / (2*DY)

alpha_c = - data['VT']['pow_c']/(data['V']['pow_c']**0.5 * data['T']['pow_c']**0.5)

lat_k = tile(s.lat[:,newaxis], (1, s.Nk))    
lat_c = tile(s.lat[:,newaxis], (1, Nc))
lat_om = s.lat

# put clim info into the array
data['V']['pow_k_clim'] = [-4.5,-2.5]
data['V']['pow_om_clim'] = [-5,-3]
data['V']['pow_c_clim'] = [-5,-3]
data['V']['log'] = True
data['V']['units'] = r'm$^2$ s$^{-2}$'
data['V']['cmap'] = get_cmap('CMRmap_r')
data['V']['title'] = r'log$_{10}$($\overline{|V|^2}$)'
data['U']['pow_k_clim'] = [-4.5,-2.5]
data['U']['pow_om_clim'] = [-5,-3]
data['U']['pow_c_clim'] = [-5,-3]
data['U']['log'] = True
data['U']['units'] = r'm$^2$ s$^{-2}$'
data['U']['cmap'] = get_cmap('CMRmap_r')
data['U']['title'] = r'log$_{10}$($\overline{|U|^2}$)'
data['T']['pow_k_clim'] = [-3,-1]
data['T']['pow_om_clim'] = [-3,-1]
data['T']['pow_c_clim'] = [-3,-1]
data['T']['log'] = True
data['T']['units'] = r'K$^2$'
data['T']['cmap'] = get_cmap('CMRmap_r')
data['T']['title'] = r'log$_{10}$($\overline{|\Theta|^2}$)'
data['VT']['pow_k_clim'] = [-0.01, 0.01]
data['VT']['pow_om_clim'] = [-0.005, 0.005]
data['VT']['pow_c_clim'] = [-0.005, 0.005]
data['VT']['log'] = False
data['VT']['units'] = r'K m s$^{-1}$'
#data['VT']['cmap'] = get_cmap('coolwarm')
data['VT']['cmap'] = get_cmap('posneg')
data['VT']['title'] = r'$\overline{V^\ast \Theta}$'
data['VU']['pow_k_clim'] = [-6e-4, 6e-4]
data['VU']['pow_om_clim'] = [-3e-4, 3e-4]
data['VU']['pow_c_clim'] = [-3e-4, 3e-4]
data['VU']['log'] = False
data['VU']['units'] = r'm$^2$ s$^{-2}$'
#data['VT']['cmap'] = get_cmap('coolwarm')
data['VU']['cmap'] = get_cmap('posneg')
data['VU']['title'] = r'$\overline{V^\ast U}$'


dk_norm = 1e3
dom_norm = 1e5
dc_norm = 100


close('all')
#for (pow_k,pow_om,pow_c) in [(V_pow_k,V_pow_om,V_pow_c),
#               (T_pow_k,T_pow_om,T_pow_c)]:
for dname, d in data.iteritems():
    if d['log']:
        pow_k = log10(d['pow_k']/s.dk / dk_norm)
        pow_om = log10(d['pow_om']/SST.dom / dom_norm)
        pow_c = log10(d['pow_c']/dc / dc_norm)
        Ucolor = 'c-'
    else:
        pow_k = d['pow_k'] / s.dk / dk_norm
        pow_om = d['pow_om'] / SST.dom / dom_norm
        pow_c = d['pow_c'] / dc / dc_norm
        Ucolor = 'm-'
    
    fig = figure()
    clf()
    subplot(131)
    pcolormesh(s.k, lat_k, pow_k, cmap=d['cmap'], rasterized=True)
    plot((4*rdat['r_dudley'])**-1 * 2 * pi, clat, 'k-', (4*rdat['r_rossby'])**-1 * 2 * pi, clat, 'k--')
    clim(d['pow_k_clim'])
    xticks(ktick, lens)
    ylim([-60,50])
    xlim([0,5e-5])
    grid()
    title(d['title'] + r"$(\kappa)$")
    xlabel(r'$2 \pi / k$ (km)')
    ylabel('lat')
    legend([r'$L_{eddy}$',r'$L_d$'], loc='upper right')
    cb=colorbar(orientation='horizontal', extendrect=True)
    cb.ax.set_title(r'%s / 10$^{-3}$ m$^{-1}$' % d['units'])
    
    subplot(132)
    pcolormesh(SST.om, lat_om, pow_om, cmap=d['cmap'], rasterized=True)
    clim(d['pow_om_clim'])
    xticks(omtick,days)
    xlim([ -(25*day/(2*pi))**-1, (60*day/(2*pi))**-1])
    ylim([-60,50])
    grid()
    title(d['title'] + r"$(\omega)$")
    xlabel(r'$2 \pi / \omega$ (days)')
    #ylabel('lat')
    cb=colorbar(orientation='horizontal', extendrect=True)
    cb.ax.set_title(r'%s / 10$^{-5}$ s$^{-1}$' % d['units'])
    
    subplot(133)
    pcolormesh(c[:,1:-1], lat_c, pow_c[:,1:-1], cmap=d['cmap'], rasterized=True)
    clim(d['pow_c_clim'])
    plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, Ucolor)
    ylim([-60,50])
    xlim([-0.5,0.2])
    grid()
    title(d['title'] + r'$(c)$')
    xlabel(r'$c$ (m/s)')
    #ylabel('lat')
    legend([r'$c_{eddy}$',r'$c_R$',r'$U_0$'], loc='upper left')
    cb=colorbar(orientation='horizontal', extendrect=True)
    cb.ax.set_title(r'%s / 0.01 m s$^{-1}$' % d['units'])
    draw()
    fig.tight_layout()
    fig.savefig('../figures/%s/integrated_spectra_%s.pdf' % (secname,dname))
    
    

d = data['VT']

#cpowlevs = 0.2*(arange(-10,10)+0.5)*d['pow_c_clim'][1]
#cpowlevs = arange(-75,76,10)/100.
#cpowticks = arange(-7,8,2)/10.
cpowlevs = (arange(-5.,5.)+0.5)/1e3
cpowticks = arange(-4,5)/1e3

figure(figsize=(6.5,4.5))    
ax1=subplot2grid((9,1), loc=(0,0), rowspan=4)
#subplot(211)
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1] / dc_norm, cpowlevs, cmap=d['cmap'], extend='both')
#pcolormesh(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1], cmap=d['cmap'], rasterized=True)
#clim(cpowlevs[r_[0,-1]])
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'm-')
ylim([10,50]); yticks(arange(10,51,10))
ylabel('lat')
grid();
legend([r'$c_{eddy}$',r'$c_R$',r'$U_0$'], loc='upper left')
title(r'$\overline{|V^\ast \Theta|}(c)$ extra tropics (K m s$^{-1}$ / 0.01 m s$^{-1}$)')

#subplot(212)
ax2=subplot2grid((9,1), loc=(4,0), rowspan=5)
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1] / dc_norm, cpowlevs, cmap=d['cmap'], extend='both')
#pcolormesh(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1], cmap=d['cmap'], rasterized=True)
#clim(cpowlevs[r_[0,-1]])
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'm-')
xlim([-0.1,0.05])
ylim([-60,-10])
xlabel(r'$c$ (m/s)')
ylabel('lat')
ax1.set_xticks(ax2.get_xticks())
ax1.set_xticklabels([])
ax1.set_xlim([-0.1,0.05])
ax2.set_xlim([-0.1,0.05])
grid();
#tight_layout()
cb=colorbar(cax=axes((0.92,0.3,0.01,0.4)),ticks=cpowticks)
savefig('../figures/%s/VT_phase_speed_spectra_extropical.pdf' % secname)

figure(figsize=(3.25,2.5))
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1] / dc_norm, cpowlevs/10, cmap=d['cmap'], extend='both')
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'm-')
ylim([-10,10])
xlim([-1,0.5])
ylabel('lat')
xlabel(r'$c$ (m/s)')
grid(); colorbar(ticks=cpowticks/10);
title(r'$\overline{|V^\ast \Theta|}(c)$ equator (K m s$^{-1}$ / 0.01 m s$^{-1}$)')
#legend([r'$c_{eddy}$',r'$c_R$',r'$U_0$'], loc='upper left')
tight_layout()
savefig('../figures/%s/VT_phase_speed_spectra_equatorial.pdf' % secname)



# comparison of zonal averages





    
