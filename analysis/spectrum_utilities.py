import numpy as np
from datetime import datetime, timedelta
from scipy.io import netcdf
import os

earth_radius = 6.353e6

regions = {
    'Kuroshio': ((150.,180.),(30., 45.)),
    'N.E. Pacific': ((190.,220.),(30., 45.)),
    'Gulf Stream': ((300.,330.),(30., 45.)),
    'Eq. Pacific': ((190.,220.),(-7.5, 7.5)),
    'S. Pacific': ((190.,220.),(-55, -40)),
    'S. Atlantic': ((310.,340.),(-55, -40)),
}

def compute_sla_power_spectrum(region_defs,
        dset = 'unfiltered',
        min_track_length=40,
        cache_count=10
        ):
    
    if dset == 'unfiltered':
        fsuffix = 'vxxc'; Nx=320; DX=7e3
    elif dset == 'filtered':
        fsuffix = 'vfec'; Nx=160; DX=14e3
    else:
        raise ValueError('dset must be either "filtered" or "unfiltered"')    
    data_base = '/Volumes/Bucket1/Data/AVISO/global/delayed-time/along-track/%s/sla/' % dset

    # size of the wavenumber grid
    Nk = Nx/2
    K = np.fft.fftfreq(Nx, DX)
    DK = K[1] - K[0]

    start_date = datetime(2002,4,24)
    switch_date = datetime(2008,10,19)
    end_date = datetime(2012,12,31)

    DT = 10 # number of days to average over             
    Nt = (end_date - start_date).days / 10 
    fulltime,t,count = 0, 0, 0
    dates = [(start_date + timedelta(days=(5+n*DT)))
             for n in range(Nt)]
    dateo = np.array([(start_date + timedelta(days=(5+n*DT))).toordinal()
             for n in range(Nt)])

    # set up dictionary for regional results
    r = dict()
    for rname in region_defs:
        r[rname] = dict(
            sla_ps = np.zeros((Nt,Nx)),
            sla_ps_tmp = np.zeros(Nx),
            sla_ps_count = np.zeros(Nt),
            sla_ps_count_tmp = 0,
            track_cache = set(),
            dateo = dateo,
            K = K    
        )


    t = 0
    still_building_cache = True
    this_date = start_date
    while this_date <= end_date:
        if this_date > switch_date:
            jason = 'j2'
        else:
            jason = 'j1'    
        data_dir = os.path.join(data_base, jason)
        y = this_date.year
        # get this day's file
        f = 'dt_global_%s_sla_%s_%s_20140106.nc' % (
            jason,fsuffix, this_date.strftime('%Y%m%d'))
        fname = os.path.join(data_dir, '%4d' % y, f)
        try:
            nc = netcdf.netcdf_file(fname)
            lonvar = nc.variables['longitude']
            latvar = nc.variables['latitude']
            slavar = nc.variables['SLA']
            fulltime += 1
            #print f
        except IOError:
            print 'No file found for %s' % this_date 
        if (np.mod(fulltime,DT)==0):
            print this_date
            # remove noise
            #sla_ps_tmp -= sla_ps[noise_idx].mean()
            for rname in r:
                print '%30s: %g tracks' % (rname, r[rname]['sla_ps_count_tmp'])
                r[rname]['sla_ps'][t] = (
                   r[rname]['sla_ps_tmp'] / 
                   r[rname]['sla_ps_count_tmp'] )
                r[rname]['sla_ps_count'][t] = r[rname]['sla_ps_count_tmp']
                r[rname]['sla_ps_tmp'] = 0.
                r[rname]['sla_ps_count_tmp'] = 0
            t += 1
            if t == cache_count:
                print 'Track Cache Built'
                for rname in r:
                    print rname
                    print r[rname]['track_cache']
                still_building_cache = False
                
        track = nc.variables['track'][:]
        # find track boundaries
        tdiff = np.where(np.diff(track)!=0)[0]
        bnd_left = np.hstack([0,tdiff+1])
        bnd_right = np.hstack([tdiff+1, len(track)])
        Ntracks = len(bnd_left)
        track_ids = track[bnd_left]
        for n in range(Ntracks):
            tid = track_ids[n]
            # if tid in my_track_ids:
            # stopped checking track id
            lon = lonvar[bnd_left[n]:bnd_right[n]]*lonvar.scale_factor
            lat = latvar[bnd_left[n]:bnd_right[n]]*latvar.scale_factor
            for rname in region_defs:
                lon_min, lon_max = (region_defs[rname][0][0],
                                    region_defs[rname][0][1])
                lat_min, lat_max = (region_defs[rname][1][0],
                                   region_defs[rname][1][1])

                do_this_track = True
                if not still_building_cache:
                    if tid not in r[rname]['track_cache']:
                        do_this_track = False
                        
                if do_this_track:
                    # this is a slow line
                    idx = ((lon>=lon_min) & (lon<=lon_max) &
                           (lat>=lat_min) & (lat<=lat_max) )                            
                    if sum(idx) > min_track_length:
                        if still_building_cache:
                            r[rname]['track_cache'].add(tid)
                        sla = slavar.scale_factor * slavar[bnd_left[n]:bnd_right[n]][idx]
                        sla -= sla.mean()
                        sla *= np.hanning(len(sla))
                        sla_fft = np.fft.fft(sla,Nx)
                        # power spectrum
                        ps = np.real(sla_fft * sla_fft.conj()) / len(sla)
                        # check parseval's theorem
                        msla2 = np.mean(sla**2)
                        np.testing.assert_almost_equal(
                            msla2, np.mean(ps), decimal=4
                        )
                        # normalize by dk
                        r[rname]['sla_ps_tmp'] += ps/DK
                        r[rname]['sla_ps_count_tmp'] += 1
        this_date += timedelta(days=1)

    #T = np.arange(Nt)*DT/365.*12. + 4    
    #return dateo, K, sla_ps, sla_ps_count
    return r
    
    