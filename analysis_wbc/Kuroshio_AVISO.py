from pylab import *
import os
from scipy.io import netcdf
import mpl_toolkits.basemap as bm
import pyresample

data_dir = '/Volumes/Bucket1/Data/AVISO/global/delayed-time/along-track/filtered/sla/j2'
years = range(2008,2013)

# the grid we used before
area_name = 'Kuroshio Orthographic Grid'
proj_id = 'ease_ortho'
area_id = proj_id
#proj4_args = '+proj=laea +lat_0=32 +lon_0=172 +a=6371228.0 +units=m'
proj4_args = '+proj=ortho +lat_0=32 +lon_0=172 +units=m'
x_size = 400
y_size = 400
area_extent = (-1600000.,-1600000.,1600000.,1600000.)
area_def = pyresample.utils.get_area_def(area_id, area_name, proj_id, proj4_args,
                                       x_size, y_size, area_extent)

#my_track_ids = [34,110,186,8,84,160,236,58,135,210,032,108,
#             227,49,125,201,23,99,175,251,73,149,255,47,123,199]

my_track_ids = [ 75,151,227, 49,125,201, 23, 99,175,251, 73,149,225, 47,123,199, 21, 97,173,249,
                234,158, 82,  6,184,108, 32,210,134, 58,236,160, 84,  8,186,110, 34,212,136]

# define boundaries for box to examine
lat_min, lat_max = 30., 45.
lon_min, lon_max = 150.,180.

min_track_length=40

# the spacing of SLA points
DX = 14e3 
# the largest possible track length
Nx = 160
# size of the wavenumber grid
Nk = Nx/2
K = fftfreq(Nx, 14e3)

track_count = dict()

Nt = 36*5             
fulltime,t,count = 0, 0, 0
DT = 10 # number of days to average over             

sla_ft = zeros((Nt,Nx))
sla_ft_tmp = zeros(Nx)

for y in years:
    for f in os.listdir(os.path.join(data_dir, '%4d' % y)):
        try:
            nc = netcdf.netcdf_file(os.path.join(data_dir, '%4d' % y, f))
            lonvar = nc.variables['longitude']
            latvar = nc.variables['latitude']
            slavar = nc.variables['SLA']
            fulltime += 1
            print f
        except IOerror:
            break
        if (mod(fulltime,DT)==0):
            print 'Resetting'
            sla_ft[t] = sla_ft_tmp / count
            sla_ft_tmp = 0
            count = 0
            t += 1
        track = nc.variables['track'][:]
        # exclude the first and last track
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
                    sla -= sla.mean()
                    sla *= hanning(len(sla))
                    sla_ft_tmp += abs(fft(sla,Nx))
                    count += 1
                    try:
                        track_count[tid].append(len(sla))
                    except KeyError:
                        track_count[tid] = [len(sla),]

T = arange(Nt)*DT/365.*12.        
idx_winter = mod(T-3,12)<2
idx_summer = mod(T-9,12)<1

mean_spectrum_summer = (sla_ft**2)[idx_summer].mean(axis=0)
mean_spectrum_winter = (sla_ft**2)[idx_winter].mean(axis=0)

figure()
loglog(K[:Nk]*1e3, mean_spectrum_summer[:Nk])
loglog(K[:Nk]*1e3, mean_spectrum_winter[:Nk])
loglog(K[:Nk][10:30]*1e3, 1e-8*(K[:Nk][10:30]*1e3)**-(11/3.), 'k--')

EKE_ft = K**2*sla_ft**2
figure()
pcolormesh(K[1:Nk]*1e3,T,ma.masked_invalid(log10(EKE_ft[:,1:Nk])))
        