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

@pytest.mark.parametrize('count', [1, 2, 3, 4, 5, 7, 8, 9, 10, 17, 100])
def test_triangular_envelope_and_headers(count):
    import math
    bank = MicrotubeBank(count, .1, .0003, .00002,
                         header_depth_m=.002, pitch_ratio=1.2,
                         additional_internal_volume_m3=1e-7)
    pitch = bank.effective_pitch_m
    assert pitch == 1.2 * bank.outer_diameter_m
    columns = math.ceil(math.sqrt(count))
    centers = [( (i % columns + .5*((i//columns)%2))*pitch,
                 (i//columns)*math.sqrt(3)/2*pitch) for i in range(count)]
    width = max(x for x,y in centers)-min(x for x,y in centers)+bank.outer_diameter_m
    height = max(y for x,y in centers)-min(y for x,y in centers)+bank.outer_diameter_m
    dims = bank.dimensions()
    assert dims['core_width_m'] == pytest.approx(width)
    assert dims['core_height_m'] == pytest.approx(height)
    assert dims['header_gas_volume_m3'] == pytest.approx(2*width*height*.002)
    assert dims['working_gas_volume_m3'] == pytest.approx(dims['tube_gas_volume_m3']+2*width*height*.002+1e-7)
    from dataclasses import asdict
    assert MicrotubeBank(**asdict(bank)).dimensions() == dims
    for diameter in (.00008, .0003, .0018):
        changed = replace(bank, inner_diameter_m=diameter)
        assert changed.effective_pitch_m > changed.outer_diameter_m


@pytest.mark.parametrize('ratio', [0, -1, 1, float('nan'), float('inf')])
def test_invalid_pitch_ratio(ratio):
    with pytest.raises(ValueError, match='pitch_ratio'):
        MicrotubeBank(10, .1, .0003, .00002, header_depth_m=.002, pitch_ratio=ratio)


def test_pitch_inputs_are_exclusive():
    with pytest.raises(ValueError, match='pitch_ratio'):
        MicrotubeBank(10, .1, .0003, .00002, .001, .002, pitch_ratio=1.2)


def test_square_envelope_is_unchanged():
    bank = MicrotubeBank(7, .1, .0003, .00002, .001, .002)
    dims = bank.dimensions()
    assert dims['core_width_m'] == .003
    assert dims['core_height_m'] == .003
    assert dims['header_gas_volume_m3'] == 2*.003*.003*.002
