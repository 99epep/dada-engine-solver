"""Verify phase removal preserves local-frame slider motion and derivatives."""
import importlib.util
import json
from pathlib import Path
import math
import numpy as np
import pytest
from dada_solver.four_bar import FourBarLoop, CouplerOutputPoint, SliderConstraint, FourBarSliderAssembly
from dada_solver.geometry import CylinderVolumeLimits

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('compact_geometry',ROOT/'examples/compact_coupler_geometry.py')
module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def test_shared_frame_matches_original_side_states():
    path=ROOT/'outputs/compact_motor_coupler_search.json'
    records=json.loads(path.read_text())['sides']
    limits=CylinderVolumeLimits(1e-6,1e-3)
    kin=module.load_compact_coupler(path,limits,limits)
    for angle in np.linspace(0,2*math.pi,31):
        states=kin.slider_states(angle)
        for index,name in enumerate(['small','large']):
            r=records[name]; raw=r['assembly']
            assembly=FourBarSliderAssembly(FourBarLoop(**raw['loop']),CouplerOutputPoint(**raw['output']),SliderConstraint(**raw['slider']))
            t=r['crank_direction']*angle+r['crank_angle_offset']
            state=assembly.evaluate(r['crank_radius']*np.array([math.cos(t),math.sin(t)]),r['crank_direction']*r['crank_radius']*np.array([-math.sin(t),math.cos(t)]))
            assert states[index].coordinate == pytest.approx(state.coordinate,abs=1e-12)
            assert states[index].coordinate_derivative == pytest.approx(state.coordinate_derivative,abs=1e-12)
        assert states[0].crank_pin == states[1].crank_pin
