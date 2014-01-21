from pylab import *
import os
import fnmatch
from scipy.io import netcdf

#dsufs = ['20100503','20110329','20110629','20111208','20120313','20120502','20120618'] # processing date

# get dates in sync with AVISO timing
#mydate = datetime.datetime(2009,1,7) # first record of 2009
mydate = datetime.datetime(2002,6,5) 


#Nt = 52*3
# the number of weeks in the NCDC data
Nt = 486

# latitude from the GHRSST 0.25 deg dataset
Ny = 460
#Nx = 120
Nx = 200
lat = -64.875 + 0.25*arange(Ny)
lon = -179.875 + 0.25*arange(Nx)

def get_netcdf_file(date):
    ddir = os.path.join(os.environ['D'], '/Volumes/BIG/Data/AVISO/global/dt/ref/msla/merged/uv')
    dsuf = 'dt_ref_global_merged_msla_uv'
    dstr = mydate.strftime('%Y%m%d')
    fstr = '%s_%s_%s_*.nc' % (dsuf,dstr,dstr)
    for file in os.listdir(ddir):
        if fnmatch.fnmatch(file, fstr):
            print file
            return netcdf.netcdf_file(os.path.join(ddir,file))
    raise IOError('No netcdf file found for date' + date)
       

# need an AVISO netcdf file
#ncf = netcdf.netcdf_file(os.path.join(ddir, '%s_%s_%s_%s.nc' % (dsuf1,dstr,dstr,dsufs[0])))
ncf = get_netcdf_file(mydate)

lat_idx = zeros(Ny, dtype('int'))
for j in arange(Ny):
    lat_idx[j] = argmin((ncf.variables['NbLatitudes'][:]-lat[j])**2)
# AVISO lons are 1/3 degree, need to interplotate
ir = r_[540:540+90]
alon = ncf.variables['NbLongitudes'][ir] - 360

output_dir = os.path.join(os.environ['D'], 'DATASTORE.RPA','projects','cospectra',
                'AVISO_dt_ref_global_merged_msla_v','timeseries_%s' % mydate.strftime('%Y%m%d'),
                'lon%6.3fto%6.3f_lat%6.3fto%6.3f' % (lon.min(),lon.max(),lat.min(),lat.max()))

try:
    print('Using output directory:')
    print(output_dir)
    os.makedirs(output_dir)
    #for j in arange(Ny):
    #    os.makedirs(os.path.join(output_dir, 'j%03d' % j))
except OSError:
    print('output directory present')

for n in arange(Nt):
    ncf = get_netcdf_file(mydate)
    # dstr = mydate.strftime('%Y%m%d')
    # ncf = False
    # for dsuf2 in dsufs:
    #     ncfname = os.path.join(ddir, '%s_%s_%s_%s.nc' % (dsuf1,dstr,dstr,dsuf2))
    #     if os.path.exists(ncfname):
    #         ncf = netcdf.netcdf_file(ncfname)
    #         break
    # if ncf is False:
    #     raise IOError('No netcdf file found')
    # 
    v = ncf.variables['Grid_0002']
    V = ma.masked_array(v[ir], v[ir]>1e18).T[lat_idx] / 100.

    print mydate
    for j in arange(Ny):
        # interpolate 1/3 to 1/4 grid
        Vi = interp(lon,alon,V[j])
        fname = os.path.join(output_dir, 'V_j%03d.bin' % j)
        f = file(fname, 'a')
        #np.save(f, SST[j].filled(0.))
        Vi.tofile(f)
        f.close()
        
    mydate += datetime.timedelta(7)