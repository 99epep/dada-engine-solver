"""Campaign-specific motor results do not define current design defaults."""
from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / 'docs'
REMOVED = ('FOUR_STAGE_K2_SEARCH.md', 'FOUR_STAGE_OPTIMIZATION.md',
           'VALVE_PLACEMENT_SCREENING.md', 'MOTOR_TEMPERATURE_AND_VALVE_CAMPAIGN.md')


def test_consumed_campaign_notes_are_removed():
    for name in REMOVED:
        assert not (DOCS / name).exists()
    for path in DOCS.rglob('*.md'):
        assert not any(name in path.read_text() for name in REMOVED), path


def test_current_design_boundaries_remain():
    physics = (DOCS / 'PHYSICS_DECISIONS.md').read_text()
    assert 'Treat placement as a design variable' in physics
    assert 'not a\nuniversal choice' in physics
    objectives = (DOCS / 'MOTOR_RESEARCH_OBJECTIVES.md').read_text()
    assert 'frozen-design sensitivity' in objectives
    assert 'best-found evidence' in objectives
    motion = (DOCS / 'MOTION_OPTIMALITY.md').read_text()
    assert 'mechanically unrealizable' in motion
    assert 'actual six-bar kinematics' in motion
