from pylab import *
import os
from scipy.io import netcdf
import pyresample
from datetime import datetime, timedelta
from scipy.signal import periodogram
from scipy.stats import chi2

# which data set to use
#dset = 'filtered'; fsuffix = 'vfec'; Nx=160; DX=14e3 
dset = 'unfiltered'; fsuffix = 'vxxc'; Nx=320; DX=7e3

data_base = '/Volumes/Bucket1/Data/AVISO/global/delayed-time/along-track/%s/sla/' % dset

my_track_ids = [ 75,151,227, 49,125,201, 23, 99,175,251, 73,149,225, 47,123,199, 21, 97,173,249,
                234,158, 82,  6,184,108, 32,210,134, 58,236,160, 84,  8,186,110, 34,212,136]

# define boundaries for box to examine
lat_min, lat_max = 30., 45.
lon_min, lon_max = 150.,180.

min_track_length=40


# size of the wavenumber grid
Nk = Nx/2
K = fftfreq(Nx, DX)
DK = K[1] - K[0]

# Xu and Fu white noise band
noise_idx = find( (K>(35e3)**-1) & (K<(25e3)**-1))

track_count = dict()

start_date = datetime(2002,4,24)
switch_date = datetime(2008,10,19)
end_date = datetime(2012,12,31)

DT = 10 # number of days to average over             
Nt = (end_date - start_date).days / 10 
fulltime,t,count = 0, 0, 0
dates = [(start_date + timedelta(days=(5+n*DT)))
         for n in range(Nt)]
dateo = array([(start_date + timedelta(days=(5+n*DT))).toordinal()
         for n in range(Nt)])

sla_ps = zeros((Nt,Nx))
# how many full length records went into the average
sla_ps_count = zeros(Nt)
sla_ps_tmp = zeros(Nx)

this_date = start_date
while this_date <= end_date:
    if this_date >= switch_date:
        jason = 'j2'
    else:
        jason = 'j1'    
    data_dir = os.path.join(data_base, jason)
    y = this_date.year
    # get this day's file
    f = 'dt_global_%s_sla_%s_%s_20140106.nc' % (jason,fsuffix, this_date.strftime('%Y%m%d'))
    fname = os.path.join(data_dir, '%4d' % y, f)
    try:
        nc = netcdf.netcdf_file(fname)
        lonvar = nc.variables['longitude']
        latvar = nc.variables['latitude']
        slavar = nc.variables['SLA']
        fulltime += 1
        print f
    except IOError:
        print 'No file found for %s' % this_date 
    if (mod(fulltime,DT)==0):
        print 'Resetting'
        # remove noise
        #sla_ps_tmp -= sla_ps[noise_idx].mean()
        sla_ps[t] = sla_ps_tmp / count
        sla_ps_count[t] = count
        sla_ps_tmp = 0
        count = 0
        t += 1
    track = nc.variables['track'][:]
    # find track boundaries
    tdiff = find(diff(track)!=0)
    bnd_left = hstack([0,tdiff+1])
    bnd_right = hstack([tdiff+1, len(track)])
    Ntracks = len(bnd_left)
    track_ids = track[bnd_left]
    for n in range(Ntracks):
        tid = track_ids[n]
        if tid in my_track_ids:
            lon = lonvar[bnd_left[n]:bnd_right[n]]*lonvar.scale_factor
            lat = latvar[bnd_left[n]:bnd_right[n]]*latvar.scale_factor
            sla = slavar.scale_factor * slavar[bnd_left[n]:bnd_right[n]]
            idx = (lon>=lon_min) & (lon<=lon_max) & (lat>=lat_min) & (lat<=lat_max)
            if sum(idx) > min_track_length:
                sla = sla[idx]
                #sla_spec = sla_wavenumber_power_spectrum(lon, lat, sla)
                if True:
                    sla -= sla.mean()
                    sla *= hanning(len(sla))
                    sla_fft = fft(sla,Nx)
                    # power spectrum
                    ps = real(sla_fft * sla_fft.conj()) / len(sla)
                    # check parseval's theorem
                    msla2 = mean(sla**2)
                    assert abs(msla2 - mean(ps))/msla2 < 1e-10
                    sla_ps_tmp += ps
                else:
                    freq,sla_ps_tmp = periodogram(sla, DX,
                        window='hanning', nfft=Nx, detrend='linear',
                        scaling='density', return_onesided=False)
                count += 1
                try:
                    track_count[tid].append(len(sla))
                except KeyError:
                    track_count[tid] = [len(sla),]
    this_date += timedelta(days=1)

# normalize the spectrum to approximate a continuous power density
sla_ps /= (K[1] - K[0])

# normalize the power spectrum to be cm^2 / cycles / km
# (instead of m^2 / cycles / m)
sla_ps *= 0.1

month = np.array([ d.month for d in dates ])
T = arange(Nt)*DT/365.*12. + 4
      
idx_winter = where(month==1)
idx_summer = where(month==7)

mean_spectrum_summer = sla_ps[idx_summer].mean(axis=0)
std_spectrum_summer = std(sla_ps[idx_summer],axis=0)
count_summer = sla_ps_count[idx_summer].sum()
mean_spectrum_winter = sla_ps[idx_winter].mean(axis=0)
std_spectrum_winter = std(sla_ps[idx_winter],axis=0)
# remove noise
summer_noise = mean_spectrum_summer[noise_idx].mean()
winter_noise = mean_spectrum_winter[noise_idx].mean()
mean_spectrum_summer_denoised = mean_spectrum_summer - summer_noise
mean_spectrum_winter_denoised = mean_spectrum_winter - winter_noise
# seasonally varying noise
noise_model = summer_noise + (winter_noise - summer_noise)*cos(2*pi/12*(T-1))

# 95% error bars
err_low = (2*count_summer)/chi2.ppf(0.05/2,2*count_summer)
err_high = (2*count_summer)/chi2.ppf(0.95/2,2*count_summer)

##### FIGURE #####
# colors
rcParams['font.size'] = 8
c1 = rcParams['axes.color_cycle'][1]
c2 = rcParams['axes.color_cycle'][2]
figure(figsize=(3.25,2.8))
ax = subplot(111)
loglog(K[:Nk]*1e3, mean_spectrum_summer_denoised[:Nk],
    '-', color=c1, linewidth=1)
loglog(K[:Nk]*1e3, mean_spectrum_winter_denoised[:Nk],
    '-', color=c2, linewidth=1)
loglog(K[:Nk]*1e3, mean_spectrum_summer[:Nk],
    '--', color=c1, linewidth=0.5, dashes=[2,2])
loglog(K[:Nk]*1e3, mean_spectrum_winter[:Nk],
    '--', color=c2, linewidth=0.5, dashes=[2,2])
loglog(ones(2)*K[noise_idx[-1]]*1e3, [1e1,1e5],
    'k:', linewidth=0.25, dashes=[2,2])
loglog(ones(2)*K[noise_idx[0]]*1e3, [1e1,1e5],
    'k:', linewidth=0.25, dashes=[2,2])
xlim([7e-4,1e-1])
ylim([1e1,1e5])
#loglog(K[:Nk][10:25]*1e3, 5e-8*(K[:Nk][10:25]*1e3)**-(11/3.), 'k--')

# slope range
sr1 = r_[6:12]
sr2 = r_[13:22]
loglog(K[sr1]*1e3, 2e3*(K[sr1]/K[sr1].mean())**-5, 'k-', linewidth=1)
ax.text(1.7e-3,5e3,r'$k^{-5}$')
loglog(K[sr2]*1e3, 1.5e2*(K[sr2]/K[sr2].mean())**-3, 'k-', linewidth=1)
ax.text(4.2e-3,1.2e2,r'$k^{-3}$')
ax.text(K[noise_idx].mean()*1e3, 3e2, 'noise band', ha='center')

legend(['Summer','Winter'], loc='lower left', frameon=False)
xlabel('inverse wavelength (cpkm)'); ylabel(r'power density (cm$^2$ / cpkm)')
title('SLA Power Spectrum - Kuroshio')
tight_layout()
savefig('figures/sla_power_spectrum_Kuroshio.pdf')

sla_ps_denoised = sla_ps - noise_model[:,newaxis]
EKE_ft = K**2*sla_ps
figure()
#pcolormesh(K[1:Nk]*1e3,dateo,ma.masked_invalid(log10(EKE_ft[:,1:Nk])))
contourf(K[1:Nk]*1e3,dateo,
    ma.masked_invalid(log10(sla_ps[:,1:Nk])),
    arange(-6,-2,0.1),
    cmap='RdBu_r')
#xlim([250**-1, 70**-1])
gca().set_xscale('log')
gca().yaxis_date()
xlabel('K (cpkm)')
title('SLA Power Spectrum : %s' % dset)
