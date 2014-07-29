from pylab import *
import pickle

# gridded
dgr = np.load('../data/SAT_50degwide_V.npz')
# along track
dat = pickle.load(open('../analysis_wbc/global_AVISO_data.pkl'))['S. Atlantic']

# latitude range
jr = (dgr['lat'] >= -55) & (dgr['lat'] < -40)
# this one needs to be normalized by dk
EKE_spec_gr = (dgr['pow_k'][jr] / dgr['dk'][jr]).mean(axis=0)
K_gr = dgr['k'][jr].mean(axis=0)/(2*pi)
dK_gr = K_gr[1] - K_gr[0]

# along-track
Nk = 160
ssh_spec = dat['sla_ps'][:,:Nk].mean(axis=0)
K_at = dat['K'][:Nk]
dK_at = K_at[1] - K_at[0]

# remove noise
noise_idx = where( (K_at>(35e3)**-1) & (K_at<(25e3)**-1))[0]
ssh_spec -= ssh_spec[noise_idx].mean()
ssh_spec = ma.masked_less_equal(ssh_spec, 0)
EKE_spec_at = (K_at*2*pi)**2 * ssh_spec
# fix normalization
EKE_spec_at *= 1e7

# total_power
tot_pow_gr = sum(ma.masked_array( dK_gr * EKE_spec_gr,
                                 (K_gr<=1e-6) | (K_gr>=1e-5)))
tot_pow_at = sum(ma.masked_array( dK_at * EKE_spec_at,
                                 (K_at<=1e-6) | (K_at>=1e-5)))

# moments
M1_gr = sum(K_gr * EKE_spec_gr) / sum(EKE_spec_gr)
M1_at = sum(K_at * EKE_spec_at) / sum(EKE_spec_at)
M2_gr = sum( (K_gr - M1_gr)**2 * EKE_spec_gr ) / sum(EKE_spec_gr)
M2_at = sum( (K_at - M1_at)**2 * EKE_spec_at ) / sum(EKE_spec_at)

gaussian_gr = 470*exp( - (K_gr - M1_gr)**2 / (2*M2_gr) )
gaussian_at = 470*exp( - (K_at - M1_at)**2 / (2*M2_at) )


close('all')
rcParams['font.size'] = 8
figure(figsize=(3.25,2.5))
loglog(K_gr*1e3, EKE_spec_gr  / 1e3,
       K_at*1e3, EKE_spec_at  / 1e3,
#       K_gr*1e3, gaussian_gr / 1e3, 'k--',
#       K_at*1e3, gaussian_at / 1e3,  'r--'
)
ylim([1e-3,1e0])
xlim([7e-4,1e-1])
xlabel('inverse wavelength (cpkm)')
ylabel('power density (m$^2$s$^{-2}$ / cpkm)')
title('EKE Power Spectra: 55S - 40S')
legend(['gridded','along track'], loc='upper right', frameon=False)
tight_layout()
show()
savefig('../figures/gridded_vs_along-track.pdf')

