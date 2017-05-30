import numpy as np
import pysynphot
from pysynphot import observation as obs
from pysynphot import spectrum
from popstar import synthetic
import pdb


def test_binning_methods():
    """
    Compare the pysynphot binflux.sum() routine on filter integrations at different resolutions.
    Shows bug in routine: integrated flux depends on the filter resolution! 

    Fix: manually integrate the binned filter function. Shows that this method performs
    much better, getting nearly the same output flux for different filter resolutions,
    as we would expect
    """
    # We'll test an integration of the vega spectrum through the WFC3-IR F127M filter
    vega = synthetic.Vega()
    filt = pysynphot.ObsBandpass('wfc3,ir,f127m')

    # Convert to ArraySpectralElement for resampling.
    filt = spectrum.ArraySpectralElement(filt.wave, filt.throughput,
                                             waveunits=filt.waveunits)
    
    # Two rebinning schemes: one coarse and the other fine
    idx = np.where(filt.throughput > 0.001)[0]
    new_wave = np.linspace(filt.wave[idx[0]], filt.wave[idx[-1]], 1500, dtype=float)
    filt_fine = filt.resample(new_wave)

    wave_bin = vega.wave
    filt_bin = synthetic.rebin_spec(filt.wave, filt.throughput, wave_bin)        
    filt_coarse = pysynphot.ArrayBandpass(wave_bin, filt_bin)

    # Do the filter integration in 2 methods: one with pysynphot binflux,
    # the other with manual integration
    vega_obs_fine = obs.Observation(vega, filt_fine, binset=filt_fine.wave, force='taper')
    vega_obs_coarse = obs.Observation(vega, filt_coarse, binset=filt_coarse.wave, force='taper')

    fine_binflux = vega_obs_fine.binflux.sum()
    coarse_binflux = vega_obs_coarse.binflux.sum()

    diff_f = np.diff(vega_obs_fine.binwave)
    diff_f = np.append(diff_f, diff_f[-1])
    fine_manual = np.sum(vega_obs_fine.binflux * diff_f)
    
    diff_c = np.diff(vega_obs_coarse.binwave)
    diff_c = np.append(diff_c, diff_c[-1])
    coarse_manual = np.sum(vega_obs_coarse.binflux * diff_c)

    print('**************************************')
    print('Integrated flux with binflux:')
    print('fine binning: {0}'.format(fine_binflux))
    print('coarse binning: {0}'.format(coarse_binflux))
    print('And with manual integration:')
    print('fine binning: {0}'.format(fine_manual))
    print('coarse binning: {0}'.format(coarse_manual))
    print('**************************************')

    pdb.set_trace()
    return
