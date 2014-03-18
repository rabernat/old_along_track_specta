from pylab import *
import os
from scipy.io import netcdf
from scipy.interpolate import interp2d
import pyresample

# this is wider than necessary to accomodate the
# curving grid at the northern end of the domain
i0,i1 = 2890,3499

# the POP model data
popstr = 'hybrid_v5_rel04_BC5_ne120_t12_pop62'
data_dir = os.path.join('/Volumes/Bucket1/Data/', popstr)
# a sample netcdf file for loading grid information
ncfname = popstr + '.pop.h.nday1.0046-04-01.nc'
nc = netcdf.netcdf_file(os.path.join(data_dir, ncfname))

# create target grid for resampling
Ny = 460
lon = arange(-180,-130,0.1) + 0.05 # 1/10 deg. in long
lat = -64.875 + 0.25*arange(Ny) # same as the SST data
lons_out,lats_out = meshgrid(lon, lat)
output_grid_def = pyresample.geometry.GridDefinition(lons=lons_out, lats=lats_out)

# the input grid for tracer variables
tlon_in = nc.variables['TLONG'][:,i0:i1].copy()
tlon_in[tlon_in>=180] -= 360
tlat_in = nc.variables['TLAT'][:,i0:i1]
input_Tgrid_def = pyresample.geometry.GridDefinition(lons = tlon_in, lats=tlat_in)

# the input grid for velocity variables
ulon_in = nc.variables['ULONG'][:,i0:i1].copy()
ulon_in[ulon_in>=180] -= 360
ulat_in = nc.variables['ULAT'][:,i0:i1]
input_Ugrid_def = pyresample.geometry.GridDefinition(lons = ulon_in, lats=ulat_in)

# the angle for velocity vector rotation
a = nc.variables['ANGLE'][:,i0:i1]

nc.close()

# pretty fast
# it would be faster to pre-computed interp arrays, but for for some reason
# that doesn't work with GridDefinition objects (only AreaDefiniton)
#SST_out = pyresample.kd_tree.resample_nearest(input_grid_def, SST, ouput_grid_def, radius_of_influence=11000)

# daily resolutionn
DT = 1
# initial date
year0, mon0, day0 = 46, 1, 1

years = arange(46,51)
nmonths = 12
monthdays = [31,28,31,30,31,30,31,31,30,31,30,31]

varnames = ['U', 'V', 'T', 'S', 'H']
ncvnames = {'U':'U1_1', 'V':'V1_1', 'T':'SST', 'S':'SSS', 'H':'HMXL_2'}
ncvscale = {'U':0.01, 'V':0.01, 'T':1., 'S':1000., 'H':0.01}
nctgrid =  {'U': False, 'V': False, 'T': True, 'S':True, 'H':True }

output_dir = dict()
for varname in varnames:
    output_dir[varname] =  os.path.join(os.environ['D'], 'DATASTORE.RPA','projects','cospectra',
                popstr + '_' + varname,
                'lon%6.3fto%6.3f_lat%6.3fto%6.3f' % (lon.min(),lon.max(),lat.min(),lat.max()),
                'timeseries_%04d%02d%02d_%1dday' % (year0, mon0, day0, DT))
    try:
        print('Using output directory:')
        print(output_dir[varname])
        os.makedirs(output_dir[varname])
    except OSError:
        print('output directory present')


n=0
for year in years:
    for nm in range(nmonths):
        ncfname = popstr + '.pop.h.nday1.%04d-%02d-01.nc' % (year,nm+1)
        nc = netcdf.netcdf_file(os.path.join(data_dir, ncfname))
        ncvars = dict()
        for vname,ncvname in ncvnames.iteritems():
            ncvars[vname] = nc.variables[ncvname]
    
        for nd in range(monthdays[nm]):
            print "Processing " + str((year,nm,nd)) 
            data = dict()
            # get the actual data
            for vname,ncvar in ncvars.iteritems():
                data[vname] = ncvscale[vname] * np.ma.masked_equal(
                    ncvar[nd,:,i0:i1], ncvar._FillValue )
            # rotate U and V correctly through angle -a
            Utrue = np.cos(a) * data['U'] + np.sin(a) * data['V']
            Vtrue = -np.sin(a) * data['U'] + np.cos(a) * data['V']
            data['U'] = Utrue
            data['V'] = Vtrue
        
            # now regrid each field
            for vname, d in data.iteritems():
                if nctgrid[vname]:
                    input_grid_def = input_Tgrid_def
                else:
                    input_grid_def = input_Ugrid_def    
                dout = pyresample.kd_tree.resample_nearest(
                    input_grid_def, d, output_grid_def, radius_of_influence=11000)
            
                # loop through j range and write output
                for j in arange(Ny):
                    fname = os.path.join(output_dir[vname], '%s_j%03d.f4.bin' % (vname,j))
                    if n==0:
                        # overwrite if it is the first entry
                        f = file(fname, 'w')
                    else:
                        # otherwise append
                        f = file(fname, 'a')
                    dout[j].filled(0.).astype('f4').tofile(f)
                    f.close()
            n += 1
