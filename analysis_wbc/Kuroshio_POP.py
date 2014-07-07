from pylab import *
import os
from scipy.io import netcdf
import mpl_toolkits.basemap as bm
import pyresample

# this is wider than necessary to accomodate the
# curving grid at the northern end of the domain
i0,i1 = 2600,3100
j0,j1 = 1325,1725

# the POP model data
popstr = 'hybrid_v5_rel04_BC5_ne120_t12_pop62'
data_dir = os.path.join('/Volumes/Bucket1/Data/', popstr)
output_dir = os.path.join(data_dir, 'processed')
# a sample netcdf file for loading grid information
ncfname = popstr + '.pop.h.nday1.0046-04-01.nc'
nc = netcdf.netcdf_file(os.path.join(data_dir, ncfname))

# the input grid for tracer variables
tlon = nc.variables['TLONG'][j0:j1,i0:i1].copy()
tlon[tlon>=180] -= 360
tlat = nc.variables['TLAT'][j0:j1,i0:i1]
input_Tgrid_def = pyresample.geometry.GridDefinition(lons = tlon, lats=tlat)

# the input grid for velocity variables
ulon = nc.variables['ULONG'][j0:j1,i0:i1].copy()
ulon[ulon>=180] -= 360
ulat = nc.variables['ULAT'][j0:j1,i0:i1]
input_Ugrid_def = pyresample.geometry.GridDefinition(lons = ulon, lats=ulat)

# the angle for velocity vector rotation
a = nc.variables['ANGLE'][j0:j1,i0:i1]

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

# set up output file for writing
f = netcdf.netcdf_file(os.path.join(output_dir, 'Kuroshio_Ext.nc'), 'w', version=2)# mmap=False)
f.history = 'Processed from hybrid_v5_rel04_BC5_ne120_t12_pop62'

f.createDimension('time', sum(monthdays)*len(years))
ftime = f.createVariable('time', 'i', ('time',))
ftime.units = nc.variables['time'].units

f.createDimension('x', x_size)
fx = f.createVariable('x', 'f', ('x',))
fx[:] = area_def.proj_x_coords
fx.units = 'm'
f.createDimension('y', y_size)
fy = f.createVariable('y', 'f', ('y',))
fy[:] = area_def.proj_y_coords
fy.units = 'm'

output_variables = dict()
for shortv,v in ncvnames.iteritems():
    myv = f.createVariable(v, 'f', ('time','y','x'))
    myv.units = nc.variables[v].units
    output_variables[shortv] = myv
    
nc.close()

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
            data_kd = dict()
            
            # get the actual data
            for vname,ncvar in ncvars.iteritems():
                data[vname] = ncvscale[vname] * np.ma.masked_equal(
                    ncvar[nd,j0:j1,i0:i1], ncvar._FillValue )
            
            # rotate U and V correctly through angle -a
            Utrue = np.cos(a) * data['U'] + np.sin(a) * data['V']
            Vtrue = -np.sin(a) * data['U'] + np.cos(a) * data['V']
            data['U'] = Utrue
            data['V'] = Vtrue
            
            # resample
            for vname in ncvars.iterkeys():
                if nctgrid[vname]:
                    source_grid = input_Tgrid_def
                else:
                    source_grid = input_Ugrid_def
                data_kd[vname] = pyresample.kd_tree.resample_gauss(input_Tgrid_def, data[vname],
                                           area_def, radius_of_influence=50000, neighbours=10,
                                           sigmas=6000, fill_value=0.)
                output_variables[vname][n] = data_kd[vname]
            #f.sync()
            n += 1
            
            
f.flush()
#f.close()

