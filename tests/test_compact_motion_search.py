"""Check synthesis sampling against the established slider backend."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('compact_search', ROOT/'examples/search_compact_motor_motion.py')
search = importlib.util.module_from_spec(spec)
spec.loader.exec_module(search)


def test_reference_slider_normalization_and_analytic_derivative():
    kin = build_model(load_simulation_configuration(ROOT/'examples/cooling_cell_x2_reference.toml')).kinematics
    angles = np.linspace(0, 2*np.pi, 1440, endpoint=False)
    q, dq, metrics = search.evaluate(kin.small_assembly, kin.crank_radius, kin.crank_direction, kin.crank_angle_offset, angles)
    numerical = (np.roll(q,-1)-np.roll(q,1))/(2*(angles[1]-angles[0]))
    assert np.max(np.abs(numerical-dq)) < 1e-4
    assert metrics['stroke_to_envelope'] > .12
    assert metrics['minimum_transmission_sine'] > .5
    unchanged, phase = search.candidate(kin.small_assembly, [1,1,1,1,1,0,0])
    assert unchanged == kin.small_assembly
    assert phase == 0


def test_coupler_output_uses_finite_rod_and_correct_derivative():
    from dada_solver.four_bar import CouplerOutputPoint
    kin = build_model(load_simulation_configuration(ROOT/'examples/cooling_cell_x2_reference.toml')).kinematics
    assembly, phase = search.candidate(kin.small_assembly, [1,1,.5,.1,1,0,0], 'coupler')
    assert isinstance(assembly.output, CouplerOutputPoint)
    assert assembly.slider == kin.small_assembly.slider
    angles = np.linspace(0, 2*np.pi, 1440, endpoint=False)
    q, dq, _ = search.evaluate(assembly, kin.crank_radius, kin.crank_direction, kin.crank_angle_offset, angles)
    numerical = (np.roll(q,-1)-np.roll(q,1))/(2*(angles[1]-angles[0]))
    assert np.max(np.abs(numerical-dq)) < 1e-4
    assert assembly.output.along_coupler == pytest.approx(.5*assembly.loop.coupler_length)
