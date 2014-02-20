from pylab import *
import os
import fnmatch 
import re

# how should this work?
# sec = Sector()
# sec.timeseries_sets
# # a dict telling us which timeseries are available
# Nc = 100.
# for j in arange(sec.Ny):
#     ts = sec.timeseries_sets['NCDC_AVHRR_AMSR_OI-1day'].load(j):
#     subplot(121)
#     ts.pcolor_hovmueller()
#     subplot(122)
#     ts.pcolor_spectrum()
#     T2f_k[j] = ts.sum_in_k()
#     T2f_om[j] = ts.sum_in_om()
#     T2f_c[j],C[j] = ts.sum_in_c(Nc)

class Sector:
    
    def __init__(self, Nx=200, Nxdata=None, Ny=460, dlon=0.25, dlat=0.25,
            lonmin=-179.875, latmin=-64.875):
        
        if Nxdata is None:
            Nxdata = Nx
        
        self.Nx = Nx
        self.Nxdata = Nxdata # the actual size of the file
        self.Ny = Ny
        # coordinates
        self.lat = latmin + dlat*arange(Ny)
        self.lon = lonmin + dlon*arange(Nx)
        self.londata = lonmin + dlon*arange(Nxdata)
        # physical dimensions
        A = 6371 * 1000. # earth radius
        self.dX = (self.lon[1] - self.lon[0])* cos(self.lat*pi/180.) * 2*pi*A / 360
        self.x = arange(-Nx/2,Nx/2)[newaxis,:] * self.dX[:,newaxis]
        self.L = self.dX * Nx
        # wavenumbers
        self.Nk = Nx/2
        self.k = zeros((Ny,self.Nk))
        for j in arange(Ny):
            self.k[j]  = 2*pi * fftshift(fftfreq(Nx,self.dX[j]))[self.Nk:]
        self.dk = diff(self.k,axis=1)[:,r_[0,0:self.Nk-1]]
            
        # where all the temporary data lives
        self.base_dir = os.path.join(os.environ['D'], 'DATASTORE.RPA','projects','cospectra')
        self.timeseries_sets = dict()
            
    def search_for_timeseries_data(self):
        base_dir = self.base_dir
        sector_prefix = 'lon%6.3fto%6.3f_lat%6.3fto%6.3f' % (
                            self.londata.min(),self.londata.max(),self.lat.min(),self.lat.max())
        for ddir in os.listdir(base_dir):
            for secdir in os.listdir(os.path.join(base_dir,ddir)):
                if secdir==sector_prefix:
                    for tsdir in os.listdir(os.path.join(base_dir,ddir,secdir)):
                        if fnmatch.fnmatch(tsdir, 'timeseries_*_*day'):
                            self.add_timeseries_set(ddir,
                                os.path.join(base_dir,ddir,secdir,tsdir))
    
    def add_timeseries_set(self, dset_name, tsdir):
        print('dset_name: %s | tsdir: %s' % (dset_name, tsdir))
        tsstr = os.path.basename(tsdir)
        crap,datestr,daystr = tsstr.split('_')
        tsname = dset_name + '-' + datestr + '_' + daystr
        dTday = int(daystr[:-3])
        year,month,day = int(datestr[:4]), int(datestr[4:6]), int(datestr[6:8])
        self.timeseries_sets[tsname] = TimeSeriesSet(sector=self,
            start_date=datetime.datetime(year,month,day), dTday=dTday, data_dir=tsdir)
        print('Added timeseries ' + tsname)

class TimeSeriesSet:
    
    def __init__(self, sector, start_date, dTday, data_dir):
        self.sector = sector
        self.start_date = start_date
        self.data_dir = data_dir
        self.dTday = dTday
        # figure out prefix
        data_fname = os.listdir(data_dir)[0]
        print data_fname
        var_str,dtype_str,suf = data_fname.split('.')
        if not dtype_str in ['f4', 'f8']:
            raise ValueError('Did not find properly formatted data files in ' + data_dir)
        self.dtype_str = dtype_str
        self.data_prefix = var_str.split('_')[0]
            
    def load(self, j, Nt=None, **kwargs):
        fname = os.path.join(self.data_dir, '%s_j%03d.%s.bin' % (self.data_prefix,j,self.dtype_str) )
        return TimeSeries(fname, dtype(self.dtype_str), self.sector, j, self.dTday, Nt, **kwargs)        
        
class TimeSeries:
    
    def __init__(self, filename, dtype, sector, j, dTday, Nt=None,
            remove_temporal_mean=True, remove_zonal_mean=False, ft_normfac=1):
        
        self.sector = sector
        self.dTday = dTday
        self.dT = dTday*24*60*60.
        self.k = self.sector.k[j]
        self.lat = self.sector.lat[j]
        self.L = self.sector.L[j]
        self.dk = self.sector.dk[j]
        self.dX = self.sector.dX[j]
        
        self.ts_data = fromfile(filename, dtype=dtype)
        Ntreal = len(self.ts_data) / self.sector.Nxdata
        self.ts_data.shape = (Ntreal,self.sector.Nxdata)
        # truncate in longitude if necessary
        self.ts_data = self.ts_data[:,:self.sector.Nx]
        
        if Nt is None:
            self.Nt = Ntreal
        else:
            # truncate
            self.Nt = Nt
            self.ts_data = self.ts_data[:self.Nt]
            
        # remember we have to define omega as negative
        self.om = -2*pi * fftshift(fftfreq(self.Nt,self.dT))
        self.per = self.dT * self.Nt
        self.t = arange(self.Nt) * self.dT
        self.dom = abs(self.om[1] - self.om[0])
        
        if remove_zonal_mean:
            Tp = self.ts_data - self.ts_data.mean(axis=1)[:,newaxis]
        else:
            # the data minus the zonal and temporal mean
            Tp = self.ts_data - self.ts_data.mean()
        
        if remove_temporal_mean:
            Tp = Tp - Tp.mean(axis=0)[newaxis,:]
        
        # this sends the filtered data back to the parent
        self.ts_data_filtered = Tp
        
        # calculate wavenumber frequency spectrum
        self.ft_data = ft_normfac * fftshift(fftn(Tp))[:,self.sector.Nk:]
        # parseval's theorem: the integral of the square in x,t = integral of the square in k,om
        #self.intTp2 = sum(Tp**2 * self.sector.dX[j] * self.dT)
        # normalize
        #self.ft_data = self.ft_data * (self.intTp2 / sum(real(self.ft_data*conj(self.ft_data)) * self.dk[:,newaxis] * self.dom))**0.5
    
    def tot_power(self):
        return self.sum_over_om(self.power_in_om())
        
    def power_in_k(self):
        """Power in wavenumber space.
        The sum of all components of the same wavenumber k."""
        return self.sum_over_om(real(self.ft_data * self.ft_data.conj()))
        
    def sum_over_om(self,field):
        return field.sum(axis=0) / self.Nt
    
    def power_in_om(self):
        """Power in frequency space.
        The sum of all components of the same frequency om."""
        return self.sum_over_k(real(self.ft_data * self.ft_data.conj()))
        
    def sum_over_k(self,field):
        # this weird stuff is necessary to account for the fact that we only
        # keep half of the spectrum. The form used here satisfies Parveval integrals
        return (0.5*field[:,0] + field[:,1:].sum(axis=1)) / self.sector.Nk

    def power_in_c(self, Nc=101):
        """Power density in phase speed space.
        The sum of all components of the same phase speed c = om/k.
        Also returns dc, the size of the c-bins. Since this is not uniform,
        care must be taken in plotting."""
        return self.sum_in_c(real(self.ft_data * self.ft_data.conj()),Nc)       
    
    def sum_in_c(self, field, Nc=101):
        # should automatically mask the zero wavenumber
        #C = ma.masked_invalid(self.om[:,newaxis] / self.k[newaxis,:])
        C = self.om[:,newaxis] / self.k[newaxis,:]
        self.C = C
        
        # minimum phase speed is a wave that will cross the sector over the period
        Cmin = self.L / self.per / 10
        Cmax = abs(ma.masked_invalid(C)).max()
        
        if not mod(Nc,2):
            raise ValueError('Nc should be an odd number')

        # old method based on log spacing (ad hoc, incorrect)
        #c = hstack( [-inf,
        #        -logspace(log2(Cmin),log2(Cmax),(Nc+1)/2,base=2)[::-1],
        #        logspace(log2(Cmin),log2(Cmax),(Nc+1)/2,base=2), inf ] )
        
        # new method based on equal angle spacing
        cstar = self.om.max() / self.k.max()
        c = -cstar/tan(linspace(self.sector.Nk**-1,pi-self.sector.Nk**-1,Nc+1))
        c = hstack( [-inf, c, inf] )

        dc = diff(c)
        T2f_c = zeros(Nc+2)
        Cpts = zeros(Nc+2)
        
        # need to account for the weird k stuff
        kmask = ones(field.shape, dtype='float')
   
        Cpts_prev = 0
        T2f_c_prev = 0
        for i in arange(Nc+1):        
            #idx = ( C >= c[i] ) & ( C < c[i+1] )
            #Cpts[i] = idx.sum()
            #T2f_c[i] = (kmask*field)[idx].sum()
            idx = (C <= c[i+1])
            Cpts_new = idx.sum()
            Cpts[i] = Cpts_new - Cpts_prev
            Cpts_prev = Cpts_new
            T2f_c_new = (kmask*field)[idx].sum()
            T2f_c[i] = T2f_c_new - T2f_c_prev
            T2f_c_prev = T2f_c_new

        # need to divide by Cpts to get proper normalization?
        return T2f_c, 0.5*(c[1:]+c[:-1]), dc, Cpts
        
        
    def convert_om_to_c(self, field, Nc=101):
        # should automatically mask the zero wavenumber
        #C = ma.masked_invalid(self.om[:,newaxis] / self.k[newaxis,:])
        C = self.om[:,newaxis] / self.k[newaxis,:]
        self.C = C
        
        # minimum phase speed is a wave that will cross the sector over the period
        Cmin = self.L / self.per / 10
        Cmax = abs(ma.masked_invalid(C)).max()
        
        # new method based on equal angle spacing
        cstar = self.om.max() / self.k.max()
        c = -cstar/tan(linspace(self.sector.Nk**-1,pi-self.sector.Nk**-1,Nc))

        field_c = zeros((Nc, self.sector.Nk))
        for i in range(self.sector.Nk):
            field_c[:,i] = interp(c, C[:,i][::-1], field[:,i][::-1],
                            left=nan, right=nan)

        return field_c, c
