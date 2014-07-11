from pylab import *
import os
from scipy.io import netcdf
import pyresample
from datetime import datetime, timedelta
from scipy.signal import periodogram

# which data set to use
dset = 'filtered'; fsuffix = 'vfec'; Nx=160; DX=14e3 
#dset = 'unfiltered'; fsuffix = 'vxxc'; Nx=320; DX=7e3

data_dir = '/Volumes/Bucket1/Data/AVISO/global/delayed-time/along-track/%s/sla/j2' % dset

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

track_count = dict()

start_date = datetime(2008,10,19)
end_date = datetime(2012,12,31)

DT = 10 # number of days to average over             
Nt = (end_date - start_date).days / 10 
fulltime,t,count = 0, 0, 0
dates = [(start_date + timedelta(days=(5+n*DT)))
         for n in range(Nt)]
dateo = array([(start_date + timedelta(days=(5+n*DT))).toordinal()
         for n in range(Nt)])

sla_ps = zeros((Nt,Nx))
sla_ps_tmp = zeros(Nx)

this_date = start_date
while this_date <= end_date:
    y = this_date.year
    # get this day's file
    f = 'dt_global_j2_sla_%s_%s_20140106.nc' % (fsuffix, this_date.strftime('%Y%m%d'))
    fname = os.path.join(data_dir, '%4d' % y, f)
    try:
        nc = netcdf.netcdf_file(fname)
        lonvar = nc.variables['longitude']
        latvar = nc.variables['latitude']
        slavar = nc.variables['SLA']
        fulltime += 1
        print f
    except IOerror:
        print 'No file found for %s' % this_date 
    if (mod(fulltime,DT)==0):
        print 'Resetting'
        sla_ps[t] = sla_ps_tmp / count
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

# normalize the power spectrum to be cm^2 / cycles / km
# (instead of m^2 / cycles / m)
sla_ps *= 0.1


T = arange(Nt)*DT/365.*12.        
idx_winter = array([ d.month==1 for d in dates ])
idx_summer = array([ d.month==7 for d in dates ])

mean_spectrum_summer = sla_ps[idx_summer].mean(axis=0)
mean_spectrum_winter = sla_ps[idx_winter].mean(axis=0)

figure()
loglog(K[:Nk]*1e3, mean_spectrum_summer[:Nk])
loglog(K[:Nk]*1e3, mean_spectrum_winter[:Nk])
#loglog(K[:Nk][10:25]*1e3, 5e-8*(K[:Nk][10:25]*1e3)**-(11/3.), 'k--')
#loglog(K[:Nk][8:20]*1e3, 1e-12*(K[:Nk][8:20]*1e3)**-5, 'k-')

legend(['Summmer','Winter',r'$K^{-11/3}$'])
xlabel('K (cpkm)'); ylabel('Power')
title('SLA Power Spectrum : %s' % dset)

EKE_ft = K**2*sla_ps
figure()
pcolormesh(K[1:Nk]*1e3,dateo,ma.masked_invalid(log10(EKE_ft[:,1:Nk])))
#xlim([250**-1, 70**-1])
#gca().set_xscale('log')
gca().yaxis_date()
xlabel('K (cpkm)')
title('SLA Power Spectrum : %s' % dset)
