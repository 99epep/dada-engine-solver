"""Explicit exchanger construction boundary, independent of channel geometry."""
from dataclasses import dataclass, replace
import math
from typing import Protocol, Mapping, runtime_checkable
from dada_solver.hydraulics import HydraulicFlowModel
from dada_solver.heat_transfer import HeatTransferModel


class LumpedWallThermalModel(Protocol):
    """Optional one-wall-energy capability of the current dynamic integrator.

    Rates are heat into gas, heat from external source, and wall energy change
    in W. The legacy air_heat_w key denotes external-source heat, not a required
    external fluid. Models with more storage states need a later integrator.
    """
    wall_capacity_j_k: float

    def rates(self, gas_temperature_k: float, wall_energy_j: float) -> Mapping[str, float | None]: ...


@dataclass(frozen=True, slots=True)
class ExchangerComponents:
    """Derived component closure; never independent optimization coordinates."""
    gas_volume_m3: float
    inlet: HydraulicFlowModel
    outlet: HydraulicFlowModel
    heat_transfer: HeatTransferModel | None = None
    wall_thermal: LumpedWallThermalModel | None = None
    metadata: tuple[tuple[str, object], ...] = ()
    validity_domain: tuple[str, ...] = ()

    def __post_init__(self):
        if not math.isfinite(self.gas_volume_m3) or self.gas_volume_m3 <= 0:
            raise ValueError('Exchanger gas volume must be finite and positive.')
        if (self.heat_transfer is None) == (self.wall_thermal is None):
            raise ValueError('Choose exactly one static or lumped-wall thermal closure.')


@runtime_checkable
class ExchangerModel(Protocol):
    """A selected family derives all thermodynamic components from its design."""
    def build(self) -> ExchangerComponents: ...


def connect_exchangers(model, heat_in: ExchangerModel, heat_out: ExchangerModel):
    """Connect independently constructed families outside the integration loop.

    The graph, two gas storage nodes and passive outlet valves are retained.
    Current integrations support two static closures or two one-wall closures.
    """
    # Placement belongs to the connected thermodynamic model.  Microtube
    # designs expose this optional field; other exchanger families remain
    # untouched and keep using their generic interface.
    if hasattr(heat_in, 'valve_placement'):
        heat_in = replace(heat_in, valve_placement=model.heat_in_valve_placement)
    if hasattr(heat_out, 'valve_placement'):
        heat_out = replace(heat_out, valve_placement=model.heat_out_valve_placement)
    incoming, outgoing = heat_in.build(), heat_out.build()
    changed = replace(model,
        machine_volumes=replace(model.machine_volumes,
            cold_heat_exchanger=incoming.gas_volume_m3,
            hot_heat_exchanger=outgoing.gas_volume_m3),
        small_cold_link=incoming.inlet, large_hot_link=outgoing.inlet,
        cold_large_valve=replace(model.cold_large_valve, flow_model=incoming.outlet),
        hot_small_valve=replace(model.hot_small_valve, flow_model=outgoing.outlet))
    if incoming.wall_thermal is not None and outgoing.wall_thermal is not None:
        from dada_solver.exchangers.air_wall import AirWallMotor
        return AirWallMotor(changed, incoming.wall_thermal, outgoing.wall_thermal)
    if incoming.heat_transfer is not None and outgoing.heat_transfer is not None:
        return replace(changed, cold_heat_transfer=incoming.heat_transfer,
                       hot_heat_transfer=outgoing.heat_transfer)
    raise ValueError('Mixed static/wall storage requires a dedicated state-layout integrator.')
