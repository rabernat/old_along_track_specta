from pylab import *
import h5py
import os
import mycolors
from scipy.ndimage.filters import gaussian_filter1d, gaussian_filter

# which data to use
secname = '50degwide'
prefix = 'SAT_%s' % secname
# the different variables available
varnames = ['V','U','T','VT','VU']
# load data
data = dict()
for v in varnames:
    data[v] = dict(np.load('../data/%s_%s.npz' % (prefix, v)))
# load grid info from data
d = data['T']
k, om, c = d['k'], d['om'], d['c']
dk, dom, dc = d['dk'], d['dom'], d['dc']
Nk, Nom, Nc = k.shape[1], len(om), c.shape[1]
lat = d['lat']
Nlat = len(lat)
lat_k = tile(lat[:,newaxis], (1, Nk))    
lat_c = tile(lat[:,newaxis], (1, Nc))
lat_om = lat

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

# figure stuff
day = 24*60*60.
days = array([-15,-30,-60,inf,60,30,15])
omtick = (day * days.astype('f4') / 2 / pi)**-1
lens =  array([1000,350,200,100,80,50])
ktick = (1000. * lens.astype('f4') / 2 / pi)**-1

rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 'small'
rcParams['axes.formatter.limits'] = [-2, 2]

# plotting configuration info
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
data['VT']['cmap'] = get_cmap('posneg')
data['VT']['title'] = r'$\overline{V^\ast \Theta}$'
data['VU']['pow_k_clim'] = [-6e-4, 6e-4]
data['VU']['pow_om_clim'] = [-3e-4, 3e-4]
data['VU']['pow_c_clim'] = [-3e-4, 3e-4]
data['VU']['log'] = False
data['VU']['units'] = r'm$^2$ s$^{-2}$'
data['VU']['cmap'] = get_cmap('posneg')
data['VU']['title'] = r'$\overline{V^\ast U}$'

dk_norm = 1e3
dom_norm = 1e5
dc_norm = 100

close('all')
for dname, d in data.iteritems():
    if d['log']:
        pow_k = log10(d['pow_k']/ dk / dk_norm)
        pow_om = log10(d['pow_om']/ dom / dom_norm)
        pow_c = log10(d['pow_c']/ dc / dc_norm)
        Ucolor = 'c-'
    else:
        pow_k = d['pow_k'] / dk / dk_norm
        pow_om = d['pow_om'] / dom / dom_norm
        pow_c = d['pow_c'] / dc / dc_norm
        Ucolor = 'm-'
    
    fig = figure(figsize=(6.5,4.5))
    clf()

    subplot(131)
    pcolormesh(k, lat_k, pow_k, cmap=d['cmap'], rasterized=True)
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
    cb.ax.set_title(r'%s / 10$^{-3}$ m$^{-1}$' % d['units'],
        {'fontsize': rcParams['axes.labelsize'],
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'})
    setp(cb.ax.get_xticklabels()[::2], visible=False)
    
    subplot(132)
    pcolormesh(om, lat_om, pow_om, cmap=d['cmap'], rasterized=True)
    clim(d['pow_om_clim'])
    xticks(omtick,days)
    xlim([ -(25*day/(2*pi))**-1, (60*day/(2*pi))**-1])
    ylim([-60,50])
    grid()
    title(d['title'] + r"$(\omega)$")
    xlabel(r'$2 \pi / \omega$ (days)')
    #ylabel('lat')
    cb=colorbar(orientation='horizontal', extendrect=True)
    cb.ax.set_title(r'%s / 10$^{-5}$ s$^{-1}$' % d['units'],
        {'fontsize': rcParams['axes.labelsize'],
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'})
    setp(cb.ax.get_xticklabels()[::2], visible=False)
    
    subplot(133)
    pcolormesh(c[:,1:-1], lat_c[:,1:-1], pow_c[:,1:-1], cmap=d['cmap'], rasterized=True)
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
    cb.ax.set_title(r'%s / 0.01 m s$^{-1}$' % d['units'],
        {'fontsize': rcParams['axes.labelsize'],
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'})
    setp(cb.ax.get_xticklabels()[::2], visible=False)

    draw()
    fig.tight_layout()
    fig.savefig('../figures/%s/integrated_spectra_%s.pdf' % (secname,dname))

# contours
d = data['VT']
cpowlevs = (arange(-5.,5.)+0.5)/1e3
cpowticks = arange(-4,5)/1e3

figure(figsize=(6.5,4.5))    
ax1=subplot2grid((9,1), loc=(0,0), rowspan=4)
contourf(c[:,1:-1], lat_c[:,1:-1], d['pow_c'][:,1:-1]/dc[:,1:-1] / dc_norm, cpowlevs, cmap=d['cmap'], extend='both')
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'm-')
ylim([10,50]); yticks(arange(10,51,10))
ylabel('lat')
grid();
legend([r'$c_{eddy}$',r'$c_R$',r'$U_0$'], loc='upper left')
title(r'$\overline{|V^\ast \Theta|}(c)$ extra tropics (K m s$^{-1}$ / 0.01 m s$^{-1}$)')

ax2=subplot2grid((9,1), loc=(4,0), rowspan=5)
contourf(c[:,1:-1], lat_c[:,1:-1], d['pow_c'][:,1:-1]/dc[:,1:-1] / dc_norm, cpowlevs, cmap=d['cmap'], extend='both')
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
cb=colorbar(cax=axes((0.92,0.3,0.01,0.4)),ticks=cpowticks)
savefig('../figures/%s/VT_phase_speed_spectra_extropical.pdf' % secname)

figure(figsize=(3.25,2.5))
contourf(c[:,1:-1], lat_c[:,1:-1], d['pow_c'][:,1:-1]/dc[:,1:-1] / dc_norm, cpowlevs/10, cmap=d['cmap'], extend='both')
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

show()
