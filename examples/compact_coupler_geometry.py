"""Reconstruct independently fitted sides in one physical crank frame."""
from dataclasses import replace
import json
import math
from pathlib import Path

from dada_solver.four_bar import (FourBarLoop, CouplerOutputPoint, SliderConstraint,
    FourBarSliderAssembly, SharedCrankFourBarVolumeKinematics)


def load_compact_coupler(path, small_limits, large_limits):
    report=json.loads(Path(path).read_text())
    records=[report['sides'][name] for name in ('small','large')]
    if records[0]['crank_radius'] != records[1]['crank_radius'] or records[0]['crank_direction'] != records[1]['crank_direction']:
        raise ValueError('Both sides must share crank radius and direction')
    assemblies=[]
    for r in records:
        if not r['dense_screen_passed']:
            raise ValueError('Candidate did not pass the dense geometric screen')
        raw=r['assembly']; angle=-r['crank_angle_offset']
        def rotate(x,y):
            return math.cos(angle)*x-math.sin(angle)*y, math.sin(angle)*x+math.cos(angle)*y
        loop=FourBarLoop(**raw['loop']); slider=SliderConstraint(**raw['slider'])
        px,py=rotate(loop.rocker_pivot_x,loop.rocker_pivot_y)
        sx,sy=rotate(slider.axis_origin_x,slider.axis_origin_y)
        assemblies.append(FourBarSliderAssembly(
            replace(loop,rocker_pivot_x=px,rocker_pivot_y=py),
            CouplerOutputPoint(**raw['output']),
            replace(slider,axis_origin_x=sx,axis_origin_y=sy,axis_angle=slider.axis_angle+angle)))
    return SharedCrankFourBarVolumeKinematics(
        crank_radius=records[0]['crank_radius'],small_assembly=assemblies[0],large_assembly=assemblies[1],
        small_volume_limits=small_limits,large_volume_limits=large_limits,
        small_volume_increases_with_coordinate=records[0]['volume_increases_with_coordinate'],
        large_volume_increases_with_coordinate=records[1]['volume_increases_with_coordinate'],
        crank_direction=records[0]['crank_direction'])
