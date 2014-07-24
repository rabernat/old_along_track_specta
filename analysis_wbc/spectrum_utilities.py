import numpy as np

earth_radius = 6.353e6

regions = {
    'Kuroshio': ((150.,180.),(30., 45.)),
    'N.E. Pacific': ((190.,220.),(30., 45.)),
    'Gulf Stream': ((300.,330.),(30., 45.)),
    'Eq. Pacific': ((190.,220.),(-7.5, 7.5)),
    'S. Pacific': ((190.,220.),(-55, -40)),
    'S. Atlantic': ((310.,340.),(-55, -40)),
}

def compute_sla_power_spectrum(coords,
        dset = 'unfiltered',
        min_track_length=40):
    
    if dset == 'unfiltered':
        fsuffix = 'vxxc'; Nx=320; DX=7e3
    elif dset == 'filtered':
        fsuffix = 'vfec'; Nx=160; DX=14e3
    else:
        raise ValueError('dset must be either "filtered" or "unfiltered"')    

    lat_min, lat_max = coords[1][0],coords[1][1]
    lon_min, lon_max = coords[0][0],coords[0][1]

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
        if this_date > switch_date:
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

    T = arange(Nt)*DT/365.*12. + 4    
    return dateo, sla_ps, sla_ps_count