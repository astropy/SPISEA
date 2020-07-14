# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This packages contains package tests.
"""

def test_evolution_merged_pisa_ekstrom():
    """Test the MergedPisaEkstromParsec() class in
    evolution.py
    """

    from popstar import evolution

    evo = evolution.MergedPisaEkstromParsec()

    iso1 = evo.isochrone(age=3e6)
    iso2 = evo.isochrone(age=5e6)

    # Check that the younger isochrone has a higher mass.
    assert iso1['mass'].max() > iso2['mass'].max()

    return
