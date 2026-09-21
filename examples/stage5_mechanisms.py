"""Existing exact mechanisms used only by the bounded Stage-5 experiments."""
from benchmark_solver_acceleration import ROOT
from dada_solver.six_bar import IndependentSixBarVolumeKinematics, load_six_bar_mechanism


def six_bar(limits):
    return IndependentSixBarVolumeKinematics(
        load_six_bar_mechanism(ROOT/'outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json',restart=0,second_branch=1),
        load_six_bar_mechanism(ROOT/'outputs/large_sixbar_stageL1.json'),
        limits.small_cylinder, limits.large_cylinder)
