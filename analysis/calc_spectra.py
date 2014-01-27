from pylab import *
import mycolors
import sector_analyzer
import os
import h5py
from scipy.ndimage.filters import gaussian_filter1d

# the wide sector
#Nx = 200; secname = '50degwide'
# the narrow sector
#Nx = 120; secname = '30degwide'
# should be the same
Nx = 120; Nxdata=200; secname='30degwide'

s = sector_analyzer.Sector(Nx=Nx)
s.search_for_timeseries_data()
dsets = ['AVISO_dt_ref_global_merged_msla_v-20020605_7day',
         'NCDC_AVHRR_AMSR_OI-20020605_7day']

Nc = 101
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
data = {'V':[],'T':[],'VT':[]}
for v in data.keys():
    data[v] = {'pow_k': ma.masked_array(zeros((s.Ny, Nk)),True),
               'pow_om':  ma.masked_array(zeros((s.Ny, Nt)),True),
               'pow_c':  ma.masked_array(zeros((s.Ny, Nc+2)),True),
               'cpts':  ma.masked_array(zeros((s.Ny, Nc+2)),True) }
za_data = {'Vp2':zeros(s.Ny),'Tp2':zeros(s.Ny), 'Tbar':zeros(s.Ny), 'VpTp':zeros(s.Ny)}

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
rcParams['axes.formatter.limits'] = [-3, 3]
def spectral_plot(SST,SSH):
    figure(figsize=(6.5,7.5))

    subplot(311)
    pcolormesh(SSH.om, SSH.k, log10(real(SSH.ft_data*conj(SSH.ft_data)).T), rasterized=True)
    clim([-9,-5])
    xticks(omtick,days)
    yticks(ktick,lens)
    ylim([0,6.5e-5]); xlim(array([-1,1])*5e-6)
    xlabel(r'$(2 \pi \omega)^{-1}$ (days)'); ylabel(r'$(2 \pi \kappa)^{-1}$ (km)')
    title(r'log$_{10}$ $|V|^2(\kappa,\omega)$ (m$^2$ s$^{-2}$)')
    legend([''])
    
    subplot(312)
    pcolormesh(SST.om, SST.k, log10(real(SST.ft_data*conj(SST.ft_data)).T), rasterized=True)
    xticks(omtick,days)
    yticks(ktick,lens)
    ylim([0,6.5e-5]); xlim(array([-1,1])*5e-6)
    clim([-7,-3])
    xlabel(r'$(2 \pi \omega)^{-1}$ (days)'); ylabel(r'$(2 \pi \kappa)^{-1}$ (km)')
    title(r'log$_{10}$ $|\Theta|^2(\kappa,\omega)$ (K$^2$)')
    colorbar()
    
    subplot(313)
    pcolormesh(SST.om, SST.k, real(SSH.ft_data*conj(SST.ft_data)).T, 
        cmap=get_cmap('bwr'), rasterized=True)
    xticks(omtick,days)
    yticks(ktick,lens)
    ylim([0,6.5e-5]); xlim(array([-1,1])*5e-6)
    clim(array([-1,1])*1e-5)
    xlabel(r'$(2 \pi \omega)^{-1}$ (days)'); ylabel(r'$(2 \pi \kappa)^{-1}$ (km)')
    title(r'Re$(V \Theta^\ast)$ (K m s$^{-1}$)')
    colorbar()
    
    tight_layout()
    
    savefig('../figures/individual_spectra/SST_SST_wavefreq_spectra_%g.pdf' % int(round(SST.lat)) )
    
#plot_js = arange(79,s.Ny,40)
plot_js = array([])

sstmask = zeros(s.Ny,bool)
sshmask = zeros(s.Ny,bool)
for j in arange(s.Ny):
    print(j)
    
    SSH = s.timeseries_sets[dsets[0]].load(j, remove_zonal_mean=True, Nt=Nt, ft_normfac=sqrt(2)/(Nt*s.Nx))
    SST = s.timeseries_sets[dsets[1]].load(j, remove_zonal_mean=True, Nt=Nt, ft_normfac=sqrt(2)/(Nt*s.Nx))
    
    if any(plot_js==j):
        spectral_plot(SST,SSH)
    
    myfields = []
    if any(SST.ts_data.sum(axis=0)==0.):
        sstmask[j] = True
    else:
        myfields.append(('T',  real(SST.ft_data*SST.ft_data.conj()) ))
    if any(SSH.ts_data.sum(axis=0)>1e10):
        sshmask[j] = True
    else:
        myfields.append(('V',  real(SSH.ft_data*SSH.ft_data.conj()) ))
    if (not sstmask[j]) and (not sshmask[j]):
        myfields.append(('VT', 2*real(SSH.ft_data*SST.ft_data.conj())))
    VTf = 2 * real( SSH.ft_data * SST.ft_data.conj() )
    for (v, field) in myfields:
        data[v]['pow_k'][j] = SSH.sum_over_om(field * mask)
        data[v]['pow_om'][j] = SSH.sum_over_k(field * mask)
        data[v]['pow_c'][j], c[j], dc[j], data[v]['cpts'][j] = SSH.sum_in_c(field * mask, Nc)
    Tbar = ma.masked_equal(SST.ts_data,0.).mean()
    # need to be careful how we define these
    # discard zero wavenumber, not zero frequency
    Vp = SSH.ts_data - SSH.ts_data.mean(axis=1)[:,newaxis]      
    Tp = SST.ts_data - SST.ts_data.mean(axis=1)[:,newaxis]
    # zonal and time mean
    za_data['Tbar'][j] = Tbar
    za_data['Vp2'][j] = mean(Vp**2)
    za_data['Tp2'][j] = mean(Tp**2)
    za_data['VpTp'][j] = mean(Vp*Tp)
    
mask = (sstmask | sshmask)
MHT = rho0*cp*s.L*interp(s.lat,mld_lat,mld) * ma.masked_array(
            za_data['VpTp'], mask) 
DY = (s.lat[1] - s.lat[0])*110e3
# 2nd order centered difference
MHT_smooth = MHT.filled(0.)
MHT_smooth[mask] = interp(s.lat[mask],s.lat[~mask],MHT[~mask])
MHT_smooth = gaussian_filter1d(MHT_smooth,1.5)
heating_rate = hstack([0, (MHT_smooth[2:] - MHT_smooth[:-2]), 0]) / s.L / (2*DY)

close('all')

# mht
figure(figsize=(6.5,5.5))
subplot(211)
plot(s.lat, MHT,'k',linewidth=2)
grid(); xlim([-50,50]); ylim(array([-1,1])*1.1e13)
xlabel('lat'); ylabel(r'MHT (W)')
title('Meridional Heat Transport')

# heating
subplot(212)
plot(s.lat, ma.masked_array(heating_rate, mask), 'k', linewidth=2)
grid(); xlim([-50,50]); ylim([-15,15])
xlabel('lat'); ylabel(r'Q (W m$^{-2}$)')
title('Heating Rate')
tight_layout()

savefig('../figures/%s/MHT_and_heating.pdf' % secname)

# output Tbar for advection/diffusion calc
Tbar_fine = tile( interp(arange(-80,80,0.1)+0.5, s.lat, za_data['Tbar'])[:,newaxis],[1,500] )
#Tbar_fine.astype(dtype('>f4')).tofile('../data/PACE_SST.bin')
dTdy = (za_data['Tbar'][2:]-za_data['Tbar'][:-2]) / (2*DY)

alpha_c = - data['VT']['pow_c']/(data['V']['pow_c']**0.5 * data['T']['pow_c']**0.5)

lat_k = tile(s.lat[:,newaxis], (1, s.Nk))    
lat_c = tile(s.lat[:,newaxis], (1, Nc))
lat_om = s.lat

# put clim info into the array
data['V']['pow_k_clim'] = [-2,1]
data['V']['pow_om_clim'] = [-0.5,2.5]
data['V']['pow_c_clim'] = [-3,0]
data['V']['log'] = True
data['V']['cmap'] = get_cmap('CMRmap_r')
data['V']['title'] = r'$\overline{|V^\ast V|}$'
data['T']['pow_k_clim'] = [0,2]
data['T']['pow_om_clim'] = [2,4]
data['T']['pow_c_clim'] = [-2,0]
data['T']['log'] = True
data['T']['cmap'] = get_cmap('CMRmap_r')
data['T']['title'] = r'$\overline{|\Theta ^\ast\Theta|}$'
data['VT']['pow_k_clim'] = array([-1,1])*5
data['VT']['pow_om_clim'] = array([-1,1])*5e2
data['VT']['pow_c_clim'] = array([-1,1])*0.4
data['VT']['log'] = False
#data['VT']['cmap'] = get_cmap('coolwarm')
data['VT']['cmap'] = get_cmap('posneg')
data['VT']['title'] = r'$\overline{|V^\ast \Theta|}$'

close('all')
#for (pow_k,pow_om,pow_c) in [(V_pow_k,V_pow_om,V_pow_c),
#               (T_pow_k,T_pow_om,T_pow_c)]:
for dname, d in data.iteritems():
    if d['log']:
        pow_k = log10(d['pow_k']/s.dk)
        pow_om = log10(d['pow_om']/SST.dom)
        pow_c = log10(d['pow_c']/dc)
    else:
        pow_k = d['pow_k'] / s.dk
        pow_om = d['pow_om'] / SST.dom
        pow_c = d['pow_c'] / dc
    
    fig = figure()
    clf()
    subplot(131)
    pcolormesh(s.k, lat_k, pow_k, cmap=d['cmap'], rasterized=True)
    plot((4*rdat['r_dudley'])**-1 * 2 * pi, clat, 'k-', (4*rdat['r_rossby'])**-1 * 2 * pi, clat, 'k--')
    clim(d['pow_k_clim'])
    xticks(ktick, lens)
    ylim([-60,50])
    xlim([0,6e-5])
    grid()
    title(d['title'] + r'$(\kappa)$')
    xlabel(r'$2 \pi / k$ (km)')
    ylabel('lat')
    legend([r'$L_{eddy}$',r'$L_d$'], loc='upper left')
    colorbar(orientation='horizontal')
    
    subplot(132)
    pcolormesh(SST.om, lat_om, pow_om, cmap=d['cmap'], rasterized=True)
    clim(d['pow_om_clim'])
    xticks(omtick,days)
    ylim([-60,50])
    grid()
    title(d['title'] + r'$(\omega)$')
    xlabel(r'$2 \pi / \omega$ (days)')
    ylabel('lat')
    colorbar(orientation='horizontal')
    

    subplot(133)
    pcolormesh(c[:,1:-1], lat_c, pow_c[:,1:-1], cmap=d['cmap'], rasterized=True)
    clim(d['pow_c_clim'])
    plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
    ylim([-60,50])
    xlim([-1,0.5])
    grid()
    title(d['title'] + r'$(c)$')
    xlabel(r'$c$ (m/s)')
    ylabel('lat')
    legend([r'$c_{eddy}$',r'$c_R$',r'$U_0$'], loc='upper left')
    colorbar(orientation='horizontal')
    
    draw()
    fig.tight_layout()
    fig.savefig('../figures/%s/integrated_spectra_%s.pdf' % (secname,dname))
    
    

d = data['VT']

#cpowlevs = 0.2*(arange(-10,10)+0.5)*d['pow_c_clim'][1]
cpowlevs = arange(-75,76,10)/100.
cpowticks = arange(-7,8,2)/10.
figure(figsize=(6.5,4.5))    
ax1=subplot2grid((9,1), loc=(0,0), rowspan=4)
#subplot(211)
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1], cpowlevs, cmap=d['cmap'], extend='both')
#pcolormesh(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1], cmap=d['cmap'], rasterized=True)
#clim(cpowlevs[r_[0,-1]])
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
ylim([10,50]); yticks(arange(10,51,10))
ylabel('lat')
grid();
legend([r'$c_{eddy}$',r'$c_R$',r'$U_0$'], loc='upper left')
title(r'$\overline{|V^\ast \Theta|}(c)$ (extra tropics)')

#subplot(212)
ax2=subplot2grid((9,1), loc=(4,0), rowspan=5)
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1], cpowlevs, cmap=d['cmap'], extend='both')
#pcolormesh(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1], cmap=d['cmap'], rasterized=True)
#clim(cpowlevs[r_[0,-1]])
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
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
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1]/dc[:,1:-1], cpowlevs/10, cmap=d['cmap'], extend='both')
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
ylim([-10,10])
xlim([-1,0.5])
ylabel('lat')
xlabel(r'$c$ (m/s)')
grid(); colorbar(ticks=cpowticks/10);
title(r'$\overline{|V^\ast \Theta|}(c)$ (equator)')
#legend([r'$c_{eddy}$',r'$c_R$',r'$U_0$'], loc='upper left')
tight_layout()
savefig('../figures/%s/VT_phase_speed_spectra_equatorial.pdf' % secname)



# comparison of zonal averages





    
