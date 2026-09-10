"""Lightweight composition boundary for independently selected model families."""
from dataclasses import dataclass
from dada_solver.configuration import SimulationConfiguration
from dada_solver.exchangers.base import ExchangerModel, connect_exchangers
from dada_solver.factory import build_model
from dada_solver.kinematics import KinematicsModel


@dataclass(frozen=True, slots=True)
class MachineDesign:
    """Reproducible inputs, not an optimizer or mechanical-loss model.

    Injected kinematics uses study angle, before the single operation reversal.
    If absent, configuration selects kinematics as before. With exchanger designs,
    their derived closures replace the legacy configuration's exchanger seeds.
    Campaigns must vary the designs, not those superseded seed quantities.
    """
    configuration: SimulationConfiguration
    heat_in: ExchangerModel | None = None
    heat_out: ExchangerModel | None = None
    kinematics: KinematicsModel | None = None

    def __post_init__(self):
        if (self.heat_in is None) != (self.heat_out is None):
            raise ValueError('Provide both exchanger designs or neither.')

    def build(self):
        model = build_model(self.configuration, kinematics=self.kinematics)
        if self.heat_in is None:
            return model
        return connect_exchangers(model, self.heat_in, self.heat_out)
