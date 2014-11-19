from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import spectrum_utilities

# helper function
def draw_screen_poly( lats, lons, m, **kwargs):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, **kwargs )
    return plt.gca().add_patch(poly)

lats = [ -30, 30, 30, -30 ]
lons = [ -50, -50, 50, 50 ]

# setup Lambert Conformal basemap.
m = Basemap(projection='robin',lon_0=180,resolution='c')
fig = plt.figure()
# draw coastlines.
m.drawcoastlines()
# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
#m.drawmapboundary(fill_color=plt.rcParams['axes.color_cycle'][2])
m.drawmapboundary(fill_color='w')
# fill continents, set lake color same as ocean color.
m.fillcontinents(color='0.5',lake_color='0.5')

n=0
patches = []
names = []
for r,coords in spectrum_utilities.regions.iteritems():
    lats = [coords[1][0], coords[1][1], coords[1][1], coords[1][0]]
    lons = [coords[0][0], coords[0][0], coords[0][1], coords[0][1]]
    patches.append(draw_screen_poly(lats, lons, m, edgecolor='none',
     facecolor=plt.rcParams['axes.color_cycle'][n], alpha=0.7))
    names.append(r)
    n += 1
plt.legend(patches, names, loc='lower left')
plt.show()

