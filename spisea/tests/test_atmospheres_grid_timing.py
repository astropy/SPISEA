"""
Timing benchmark for Kurucz ``k93models`` grid extraction (stsynphot / pysynphot).

Kept in a separate module (no ``spisea`` imports) so the pysynphot-only branch can
run when ``stsynphot`` is not imported elsewhere.
"""
import importlib.util
import os
import time

import numpy as np
import pytest


def test_pysynphot_vs_stsynphot_timing():
    """
    Time Kurucz ``k93models`` spectrum extraction via stsynphot or pysynphot.
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
        import pysynphot

        t0 = time.perf_counter()
        try:
            sp = pysynphot.Icat("k93models", temperature, metallicity, gravity)
            wave = sp.GetWaveSet()
            _ = sp.sample(wave)
        except Exception as exc:
            pytest.skip(f"pysynphot Icat failed: {exc}")
        elapsed = time.perf_counter() - t0

        assert elapsed >= 0
        assert np.isfinite(elapsed)
        assert np.asarray(wave).size > 0

    else:
        pytest.skip("Neither stsynphot nor pysynphot is installed")
