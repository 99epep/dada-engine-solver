"""Keep the capacity reference focused on the current transformation."""
from pathlib import Path

from dada_solver.research.rescale import extensive

DOC = Path(__file__).resolve().parents[1] / 'docs/DADA_ENGINE_RESEARCH_CAPACITY.md'


def test_capacity_role_and_current_contract():
    text = DOC.read_text()
    assert text.startswith('# Capacity scaling')
    for token in ('V2/V3', 'local_regions_v1', 'strictly positive',
                  'Constraints remain unchanged', 'ties to even', "N' L' = s N L", 'fixed_global_bounds',
                  'not geometric similarity', 'Factor `1`', 'new scientific identity'):
        assert token in text
    for name in ('DADA_ENGINE_RESEARCH.md', 'DADA_ENGINE_RESEARCH_REFERENCE.md',
                 'DADA_ENGINE_RESEARCH_COCKPIT.md', 'MICROTUBE_GAS_MODEL.md'):
        assert f']({name})' in text


def test_production_extensive_coordinate_boundary():
    assert extensive('volume.total_swept_m3')
    assert not extensive('microtube.heat_in.tube_count')
    assert not extensive('microtube.heat_in.tube_length_m')
    assert extensive('microtube.heat_in.additional_internal_volume_m3')
    assert extensive('external_stream.heat_in.mass_flow_kg_s')
    assert not extensive('operation.frequency_hz')
    assert not extensive('microtube.heat_in.pitch_ratio')


def test_no_campaign_or_cockpit_walkthrough():
    text = DOC.read_text()
    for token in ('../outputs/', 'examples/', 'human_cell_stage0',
                  '## Valve chronology availability', '## Offline table and volume plots',
                  'topology_display', '--plots volumes', 'Bounded Human Cell validation'):
        assert token not in text
