from pylab import *
import mycolors
import sector_analyzer

s = sector_analyzer.Sector()
s.search_for_timeseries_data()
dsets = ['AVISO_dt_ref_global_merged_msla_v-20020605_7day',
         'NCDC_AVHRR_AMSR_OI-20020605_7day']

Nc = 101
Nt = 486
Nk = s.Nk

# initialize results matrix
# fourier transform data
data = {'V':[],'T':[],'VT':[]}
for v in data.keys():
    data[v] = {'pow_k': zeros((s.Ny, Nk)),
               'pow_om': zeros((s.Ny, Nt)),
               'pow_c': zeros((s.Ny, Nc+2)),
               'cpts': zeros((s.Ny, Nc+2)) }
za_data = {'Vp2':zeros(s.Ny),'Tp2':zeros(s.Ny), 'Tbar':zeros(s.Ny), 'VpTp':zeros(s.Ny)}

c = zeros((s.Ny, Nc+2))
dc = zeros((s.Ny, Nc+2))

mask = ones(s.Nk)
mask[0] = 0

# figure stuff
day = 24*60*60.
days = array([-15,-30,-60,-120,120,60,30,15])
lens =  array([1000,500,200,100,50])
omtick = 2*pi*(day * days.astype('f4'))**-1
ktick = 2*pi*(1000. * lens.astype('f4'))**-1


for j in arange(s.Ny):
    print(j)
    SSH = s.timeseries_sets[dsets[0]].load(j, remove_zonal_mean=True)
    SST = s.timeseries_sets[dsets[1]].load(j, remove_zonal_mean=True)
    VTf = 2 * real( SSH.ft_data * SST.ft_data.conj() )
    for (v, field) in [ ('V',  real(SSH.ft_data*SSH.ft_data.conj()) ),
                     ('T',  real(SST.ft_data*SST.ft_data.conj()) ),
                     ('VT', 2*real(SSH.ft_data*SST.ft_data.conj())) ]:
        data[v]['pow_k'][j] = SSH.sum_over_om(field * mask)
        data[v]['pow_om'][j] = SSH.sum_over_k(field * mask)
        data[v]['pow_c'][j], c[j], dc[j], data[v]['cpts'][j] = SSH.sum_in_c(field * mask, Nc)
    Vp = SSH.ts_data - SSH.ts_data.mean()      
    Tbar = SST.ts_data.mean()
    Tp = SST.ts_data - Tbar
    # zonal and time mean
    za_data['Tbar'][j] = Tbar
    za_data['Vp2'][j] = mean(Vp**2)
    za_data['Tp2'][j] = mean(Tp**2)
    za_data['VpTp'][j] = mean(Vp*Tp)

alpha_c = - data['VT']['pow_c']/(data['V']['pow_c']**0.5 * data['T']['pow_c']**0.5)

lat_k = tile(s.lat[:,newaxis], (1, s.Nk))    
lat_c = tile(s.lat[:,newaxis], (1, Nc))
lat_om = s.lat

cdat = np.load('../aviso_mixing/data/andreas/c.npz')
clat = linspace(-80,80,160)
Udat = np.load('../aviso_mixing/data/andreas/Umean_ECCO_patch.npz')
rdat = np.load('../aviso_mixing/data/andreas/r.npz')

# put clim info into the array
data['V']['pow_k_clim'] = [2,4.2]
data['V']['pow_om_clim'] = [2,4.2]
data['V']['pow_c_clim'] = [4,6.5]
data['V']['log'] = True
data['V']['cmap'] = get_cmap('CMRmap_r')
data['V']['title'] = r'$v v^\ast$'
data['T']['pow_k_clim'] = [3,6]
data['T']['pow_om_clim'] = [4,6]
data['T']['pow_c_clim'] = [6.5,7.5]
data['T']['log'] = True
data['T']['cmap'] = get_cmap('CMRmap_r')
data['T']['title'] = r'$\theta \theta^\ast$'
data['VT']['pow_k_clim'] = array([-1,1])*2e4
data['VT']['pow_om_clim'] = array([-1,1])*2e4
data['VT']['pow_c_clim'] = array([-1,1])*3e6
data['VT']['log'] = False
#data['VT']['cmap'] = get_cmap('coolwarm')
data['VT']['cmap'] = get_cmap('posneg')
data['VT']['title'] = r'$v \theta^\ast$'

close('all')
#for (pow_k,pow_om,pow_c) in [(V_pow_k,V_pow_om,V_pow_c),
#               (T_pow_k,T_pow_om,T_pow_c)]:
for dname, d in data.iteritems():
    if d['log']:
        pow_k = log10(d['pow_k'])
        pow_om = log10(d['pow_om'])
        pow_c = log10(d['pow_c'])
    else:
        pow_k = d['pow_k']
        pow_om = d['pow_om']
        pow_c = d['pow_c']
    
    figure()
    clf()
    subplot(131)
    pcolormesh(s.k, lat_k, pow_k, cmap=d['cmap'])
    plot((4*rdat['r_dudley'])**-1 * 2 * pi, clat, 'k-', (4*rdat['r_rossby'])**-1 * 2 * pi, clat, 'k--')
    clim(d['pow_k_clim'])
    ylim([-60,50])
    xlim([0,5e-5])
    xticks(ktick, lens)
    grid()
    title(d['title'])
    xlabel(r'$2 \pi / k$ (km)')
    ylabel('lat')
    
    subplot(132)
    pcolormesh(SST.om, lat_om, pow_om, cmap=d['cmap'])
    clim(d['pow_om_clim'])
    ylim([-60,50])
    xticks(omtick,days)
    grid()
    title(d['title'])
    xlabel(r'$2 \pi / \omega$ (days)')
    ylabel('lat')

    subplot(133)
    pcolormesh(c[:,1:-1], lat_c, pow_c[:,1:-1], cmap=d['cmap'])
    clim(d['pow_c_clim'])
    plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
    ylim([-60,50])
    xlim([-1,0.5])
    grid()
    title(d['title'])
    xlabel(r'$c$ (m/s)')
    ylabel('lat')

d = data['VT']

cpowlevs = 0.1*(arange(-10,10)+0.5)*d['pow_c_clim'][1]
figure()    
subplot2grid((11,1), loc=(0,0), rowspan=4)
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1], cpowlevs, cmap=d['cmap'], extend='both')
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
xlim([-0.1,0.05])
ylim([10,50])
ylabel('lat')
grid(); colorbar()
subplot2grid((11,1), loc=(4,0), rowspan=2)
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1], cpowlevs, cmap=d['cmap'], extend='both')
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
ylim([-10,10])
xlim([-1,0.5])
ylabel('lat')
grid(); colorbar()
subplot2grid((11,1), loc=(6,0), rowspan=5)
contourf(c[:,1:-1], lat_c, d['pow_c'][:,1:-1], cpowlevs, cmap=d['cmap'], extend='both')
plot(-cdat['c_dudley'], clat, 'k-', cdat['c_doppler'], clat, 'k--', Udat['Umean_ECCO_patch'], clat, 'k:')
xlim([-0.1,0.05])
ylim([-60,-10])
xlabel(r'$c$ (m/s)')
ylabel('lat')

grid(); colorbar()

# comparison of zonal averages





    
