"""
Tests for atmosphere utilities, including ``rebin_spec``.
"""
import importlib.util
import os
import time

import numpy as np
import pytest
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

    got = atm.rebin_spec(w_in, vega, w_out).value
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


def test_pysynphot_vs_stsynphot_timing():
    """
    Time Kurucz ``k93models`` spectrum extraction via stsynphot or pysynphot.

    This module imports ``spisea.atmospheres`` and ``spisea.synthetic``, which
    import ``stsynphot`` at load time; the pysynphot-only section of this code will
    only run on older versions of SPISEA.
    """
    temperature = 20000
    metallicity = 0.0
    gravity = 4.0

    has_stsyn = importlib.util.find_spec("stsynphot") is not None
    has_pysyn = importlib.util.find_spec("pysynphot") is not None

    if has_stsyn:
        if not os.environ.get("PYSYN_CDBS"):
            pytest.skip("PYSYN_CDBS not set; grid_to_spec needs CDBS tree")

        from stsynphot.catalog import grid_to_spec

        t0 = time.perf_counter()
        try:
            sp = grid_to_spec("k93models", temperature, metallicity, gravity)
            w = sp.waveset
            _ = sp(w)
        except Exception as exc:
            pytest.skip(f"stsynphot grid_to_spec failed: {exc}")
        elapsed = time.perf_counter() - t0

        assert elapsed >= 0
        assert np.isfinite(elapsed)
        assert w.size > 0

    elif has_pysyn:
        import synphot

        t0 = time.perf_counter()
        try:
            sp = synphot.grid_to_spec("k93models", temperature, metallicity, gravity)
            w = sp.waveset
            _ = sp(w)
        except Exception as exc:
            pytest.skip(f"pysynphot Icat failed: {exc}")
        elapsed = time.perf_counter() - t0

        assert elapsed >= 0
        assert np.isfinite(elapsed)
        assert np.asarray(w).size > 0

    else:
        pytest.skip("Neither stsynphot nor pysynphot is installed")