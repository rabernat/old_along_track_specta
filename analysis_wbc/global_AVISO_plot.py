from pylab import *
import os
from scipy.stats import chi2
from datetime import datetime, timedelta
import spectrum_utilities
import pickle
from scipy.ndimage import gaussian_filter1d

result = pickle.load(open('global_AVISO_data.pkl'))

close('all')

for rname,data in result.iteritems():
    short_name = rname.replace('. ',
                    '_').replace('.','').replace(' ','_').lower()
    sla_ps = 0.1 * data['sla_ps'] # make into cm^2 / cpkm
    K = data['K']
    dateo = data['dateo']
    
    Nt, Nx = sla_ps.shape
    Nk = Nx/2
    dates = [ datetime.fromordinal(do) for do in dateo]
    month = np.array([ d.month for d in dates ])
    #T = arange(Nt)*DT/365.*12. + 4
      
    # smooth by one point in K
    sla_ps[1:-1] = (sla_ps[:-2] + sla_ps[1:-1] + sla_ps[:-2])/3.
      
    idx_winter = where(month==1)
    idx_summer = where(month==7)

    mean_spectrum_summer = sla_ps[idx_summer].mean(axis=0)
    std_spectrum_summer = std(sla_ps[idx_summer],axis=0)
    count_summer = data['sla_ps_count'][idx_summer].sum()
    mean_spectrum_winter = sla_ps[idx_winter].mean(axis=0)
    std_spectrum_winter = std(sla_ps[idx_winter],axis=0)

    # remove noise (Xu and Fu white noise band)
    noise_idx = where( (K>(35e3)**-1) & (K<(25e3)**-1))[0]
    summer_noise = mean_spectrum_summer[noise_idx].mean()
    winter_noise = mean_spectrum_winter[noise_idx].mean()
    mean_spectrum_summer_denoised = mean_spectrum_summer - summer_noise
    mean_spectrum_winter_denoised = mean_spectrum_winter - winter_noise

    # seasonally varying noise
    #noise_model = summer_noise + (winter_noise - summer_noise)*cos(2*pi/12*(T-1))

    # 95% error bars
    err_low = (2*count_summer)/chi2.ppf(0.05/2,2*count_summer)
    err_high = (2*count_summer)/chi2.ppf(0.95/2,2*count_summer)

    ##### FIGURE #####
    # colors
    rcParams['font.size'] = 8
    c1 = rcParams['axes.color_cycle'][1]
    c2 = rcParams['axes.color_cycle'][2]
    figure(figsize=(3.25,2.8))
    ax = subplot(111)
    loglog(K[:Nk]*1e3, mean_spectrum_summer_denoised[:Nk],
        '-', color=c1, linewidth=1)
    loglog(K[:Nk]*1e3, mean_spectrum_winter_denoised[:Nk],
        '-', color=c2, linewidth=1)
    loglog(K[:Nk]*1e3, mean_spectrum_summer[:Nk],
        '--', color=c1, linewidth=0.5, dashes=[2,2])
    loglog(K[:Nk]*1e3, mean_spectrum_winter[:Nk],
        '--', color=c2, linewidth=0.5, dashes=[2,2])
    loglog(ones(2)*K[noise_idx[-1]]*1e3, [1e1,1e5],
        'k:', linewidth=0.25, dashes=[2,2])
    loglog(ones(2)*K[noise_idx[0]]*1e3, [1e1,1e5],
        'k:', linewidth=0.25, dashes=[2,2])
    xlim([7e-4,1e-1])
    ylim([1e1,1e5])
    #loglog(K[:Nk][10:25]*1e3, 5e-8*(K[:Nk][10:25]*1e3)**-(11/3.), 'k--')

    if (short_name in ['kuroshio', 'gulf_stream']):
        # slope range
        sr1 = r_[6:12]
        sr2 = r_[13:22]
        loglog(K[sr1]*1e3, 2e3*(K[sr1]/K[sr1].mean())**-5, 'k-', linewidth=1)
        ax.text(1.7e-3,5e3,r'$k^{-5}$')
        loglog(K[sr2]*1e3, 1.5e2*(K[sr2]/K[sr2].mean())**-3, 'k-', linewidth=1)
        ax.text(4.2e-3,1.2e2,r'$k^{-3}$')
        ax.text(K[noise_idx].mean()*1e3, 3e2, 'noise band', ha='center')

    legend(['January','July'], loc='lower left', frameon=False)
    xlabel('inverse wavelength (cpkm)'); ylabel(r'power density (cm$^2$ / cpkm)')
    title('SLA Power Spectrum - %s' % rname)
    tight_layout()
    savefig('figures/sla_ps_seasonal_%s.pdf' % short_name)

    # sla_ps_denoised = sla_ps - noise_model[:,newaxis]
    figure(figsize=(3.25,5))
    sla_ps_denoised = sla_ps - sla_ps[:,noise_idx].mean(axis=1)[:,newaxis]
    # EKE_ft = K**2*sla_ps
    # figure()
    # #pcolormesh(K[1:Nk]*1e3,dateo,ma.masked_invalid(log10(EKE_ft[:,1:Nk])))
    ax = subplot2grid((20,1),(0,0),rowspan=16)
    axcb = subplot2grid((20,1),(19,0))
    cf=ax.contourf(K[1:Nk]*1e3,dateo,
        ma.masked_invalid(log10(sla_ps_denoised[:,1:Nk])),
        arange(1.5,5.1,0.25),
        cmap='CMRmap_r', extend='both', rasterized=True)
    ax.set_xlim([7e-4,1e-1])
    #ylim([1e1,1e5])
    # xlim([250**-1, 70**-1])
    ax.set_xscale('log')
    ax.yaxis_date()
    ax.grid()
    colorbar(cf, cax=axcb, orientation='horizontal')
    ax.set_xlabel('K (cpkm)')
    ax.set_title('SLA Spectrogram - %s' % rname)
    axcb.set_title(r'log$_{10}$ power density (cm$^2$ / cpkm)')
    subplots_adjust(bottom=0.04, hspace=0.75)
    savefig('figures/sla_ps_full_%s.pdf' % short_name)
    #tight_layout()
    
    bigfig=figure(figsize=(6.5,2.8))
    axf = bigfig.add_subplot(111)
    K0 = 25
    axf.semilogy(dateo,
        gaussian_filter1d(
            sla_ps_denoised[:,K0:K0+5].mean(axis=1),1))    
    axf.xaxis_date()
    axf.grid()
    axf.set_ylabel(r'power density (cm$^2$ / cpkm)')
    axf.set_title('SLA Spectral Power at 100 km - %s' % rname)
    savefig('figures/sla_ps_100km_timeseries_%s.pdf' % short_name)