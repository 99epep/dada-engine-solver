from dataclasses import replace

import pytest

from dada_solver.exchangers.microtube_geometry import MicrotubeBank


def test_more_shorter_tubes_preserve_core_volume_but_increase_headers():
    original = MicrotubeBank(100, .1, .00033, .0001524, .00125, .002)
    changed = replace(original, tube_count=400, tube_length_m=.025)
    a, b = original.dimensions(), changed.dimensions()
    assert a['tube_internal_area_m2'] == pytest.approx(b['tube_internal_area_m2'])
    assert a['tube_gas_volume_m3'] == pytest.approx(b['tube_gas_volume_m3'])
    assert b['header_gas_volume_m3'] == pytest.approx(4*a['header_gas_volume_m3'])
    properties = dict(density_kg_m3=3, viscosity_pa_s=2e-5)
    dp = original.laminar_tube_loss(.0005, **properties)['signed_tube_pressure_drop_pa']
    assert changed.laminar_tube_loss(.0005, **properties)['signed_tube_pressure_drop_pa'] == pytest.approx(dp/16)
    assert original.laminar_tube_loss(-.0005, **properties)['signed_tube_pressure_drop_pa'] == -dp
    assert original.laminar_tube_loss(0, **properties)['signed_tube_pressure_drop_pa'] == 0


@pytest.mark.parametrize('change', [dict(tube_count=True), dict(pitch_m=.0004),
    dict(header_depth_m=0), dict(additional_internal_volume_m3=float('nan'))])
def test_invalid_geometry(change):
    with pytest.raises(ValueError):
        replace(MicrotubeBank(100, .1, .00033, .0001524, .00125, .002), **change)
