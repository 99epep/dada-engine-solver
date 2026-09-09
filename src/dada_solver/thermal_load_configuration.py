"""TOML configuration for external cooling tasks and human power."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib

from dada_solver.thermal_loads import CoolingTask, HumanPowerEnvelope, LumpedPhaseChangeLoad


@dataclass(frozen=True, slots=True)
class CoolingCellScenario:
    load: LumpedPhaseChangeLoad
    reference_task: CoolingTask
    ambitious_task: CoolingTask
    human_power: HumanPowerEnvelope
    ambient_temperature: float


def load_cooling_cell_scenario(path: str | Path) -> CoolingCellScenario:
    """Load an external ideal cooling-cell task using SI values."""

    with Path(path).open("rb") as stream:
        data = tomllib.load(stream)
    try:
        load_data = data["thermal_load"]
        target_data = data["targets"]
        human_data = data["human_power"]
        environment = data["environment"]
        load = LumpedPhaseChangeLoad(
            material_name=str(load_data["material_name"]),
            mass=float(load_data["mass"]),
            initial_temperature=float(load_data["initial_temperature"]),
            target_temperature=float(load_data["target_temperature"]),
            initial_phase_specific_heat=float(
                load_data["initial_phase_specific_heat"]
            ),
            phase_change_specific_energy=float(
                load_data["phase_change_specific_energy"]
            ),
            phase_change_fraction=float(load_data["phase_change_fraction"]),
        )
        return CoolingCellScenario(
            load=load,
            reference_task=CoolingTask(load, float(target_data["reference_time"])),
            ambitious_task=CoolingTask(load, float(target_data["ambitious_time"])),
            human_power=HumanPowerEnvelope(
                minimum_power=float(human_data["minimum_power"]),
                nominal_power=float(human_data["nominal_power"]),
                maximum_power=float(human_data["maximum_power"]),
            ),
            ambient_temperature=float(environment["ambient_temperature"]),
        )
    except KeyError as error:
        raise ValueError(
            f"Missing cooling-cell configuration field: {error.args[0]}"
        ) from error

