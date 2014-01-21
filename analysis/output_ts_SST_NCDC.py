from pylab import *
import os
from scipy.io import netcdf


ddir = os.path.join(os.environ['D'], 'ghrsst/data/L4/GLOB/NCDC/AVHRR_AMSR_OI')
dsuf = 'NCDC-L4LRblend-GLOB-v01-fv02_0-AVHRR_AMSR_OI.nc'
year = 2009
day = 1

DT = 1
# AVISO frequency
#DT = 7 

# get dates in sync with AVISO timing
mydate = datetime.datetime(2002,6,5) 
lastdate = datetime.datetime(2011,11,4)

# weekly resolution
#Nt = 52*10
#Nt = 486 # that's all we can get
# daily resolution
Nt = floor((lastdate - mydate).days / DT) - 1

# range for analysis
#ir = r_[0:120] # 180 - 150 W
ir = r_[0:200] # 180 - 130 W
jr = r_[100:560]
Nx,Ny = len(ir),len(jr)
# latitude from the GHRSST 0.25 deg dataset
lat = -64.875 + 0.25*arange(Ny)
lon = -179.875 + 0.25*arange(Nx)

output_dir = os.path.join(os.environ['D'], 'DATASTORE.RPA','projects','cospectra',
                'NCDC_AVHRR_AMSR_OI',
                'lon%6.3fto%6.3f_lat%6.3fto%6.3f' % (lon.min(),lon.max(),lat.min(),lat.max()),
                'timeseries_%s_%1dday' % (mydate.strftime('%Y%m%d'), DT))

try:
    print('Using output directory:')
    print(output_dir)
    os.makedirs(output_dir)
    #for j in arange(Ny):
    #    os.makedirs(os.path.join(output_dir, 'j%03d' % j))
except OSError:
    print('output directory present')

for n in arange(Nt):
    year = mydate.year
    yearday = (mydate - datetime.datetime(year,1,1)).days + 1
    ncf = netcdf.netcdf_file(os.path.join(ddir, '%4d' % year, '%03d' % yearday,
                '%s-%s' % (mydate.strftime('%Y%m%d'),dsuf)))
    v = ncf.variables['analysed_sst']
    SST = ma.masked_array(v[0,jr][:,ir]*v.scale_factor + v.add_offset, (ncf.variables['mask'][0,jr][:,ir]!=1))

    print mydate
    for j in arange(Ny):
        fname = os.path.join(output_dir, 'SST_j%03d.bin' % j)
        if n==0:
            f = file(fname, 'w')
        else:
            f = file(fname, 'a')
        #np.save(f, SST[j].filled(0.))
        SST[j].filled(0.).tofile(f)
        f.close()
        
    mydate += datetime.timedelta(DT)