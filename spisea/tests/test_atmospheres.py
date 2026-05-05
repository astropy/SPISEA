"""
Tests for atmosphere utilities, including ``rebin_spec``.
"""
import numpy as np
from astropy import units as u
from synphot import Observation
from synphot import units as su
from synphot.models import Empirical1D
from synphot.spectrum import SpectralElement

from spisea import atmospheres as atm
from spisea import synthetic as syn


def _unit_flat_bandpass(waveset):
    """Unit throughput on *waveset* (same construction as ``synphot_bridge.rebin_spec``)."""
    n = len(np.asarray(waveset.to(u.AA).value))
    return SpectralElement(
        Empirical1D,
        points=waveset,
        lookup_table=np.ones(n) * su.THROUGHPUT,
    )


def test_rebin_spec_matches_synphot_observation():
    """``rebin_spec`` must return the same binned flux as ``Observation.binflux``."""
    vega = syn.vega
    w_in = vega.waveset
    w_out = w_in[::28]
    filt = _unit_flat_bandpass(w_in)

    got = atm.rebin_spec(w_in, vega, w_out)
    obs = Observation(vega, filt, binset=w_out, force="taper")
    expected = np.asarray(obs.binflux.value, dtype=np.float64).ravel()

    np.testing.assert_allclose(got, expected, rtol=0, atol=0)

    return


def test_vega_rebin_conserves_integrated_flux():
    """
    Total integrated flux (spec × unit band) is independent of binning grid.

    Uses the package Vega spectrum and the same flat bandpass as ``rebin_spec``;
    :meth:`~synphot.spectrum.BaseSourceSpectrum.integrate` must agree for fine
    and coarse ``binset`` values.
    """
    vega = syn.vega
    w_in = vega.waveset
    filt = _unit_flat_bandpass(w_in)

    obs_fine = Observation(vega, filt, binset=w_in[::3], force="taper")
    obs_coarse = Observation(vega, filt, binset=w_in[::45], force="taper")

    int_fine = obs_fine.integrate()
    int_coarse = obs_coarse.integrate()

    assert int_fine.unit == int_coarse.unit
    np.testing.assert_allclose(
        int_fine.value,
        int_coarse.value,
        rtol=1e-12,
        atol=0,
        err_msg="Integrated flux changed when only the rebin grid changed",
    )

    return


def test_rebin_spec_output_length_matches_binset():
    vega = syn.vega
    w_out = vega.waveset[::33]
    bf = atm.rebin_spec(vega.waveset, vega, w_out)
    assert bf.shape == (w_out.size,)

    return
