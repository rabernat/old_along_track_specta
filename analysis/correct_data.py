from pylab import *

prefix = 'POP_50degwide'
# the different variables available
varnames = ['V','U','T','VT','VU']
# load data
for v in varnames:
    data = dict(np.load('../data/POP_old/%s_%s.npz' % (prefix, v)))
    # renormalize data
    Nlat,Nk = data['pow_k'].shape
    Nlat,Nom = data['pow_om'].shape
    data['pow_om'] *= Nk
    data['pow_k'] *= Nom
    np.savez('../data/%s_%s.npz' % (prefix, v), **data)
    