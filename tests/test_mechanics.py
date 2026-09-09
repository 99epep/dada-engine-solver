import numpy as np
import pytest

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.mechanics import (
    PistonExternalPressures,
    PistonFaceAreas,
    extract_thermodynamic_loads,
)
from tests.test_periodic import create_static_cycle


def test_thermodynamic_mechanical_boundary_exposes_only_gas_side_loads(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, _, cycle = create_static_cycle(ideal_gas)
    areas = PistonFaceAreas(small=2.0e-3, large=4.0e-3)

    loads = extract_thermodynamic_loads(cycle, model, areas)

    np.testing.assert_allclose(loads.small_pressure, 2.0e5)
    np.testing.assert_allclose(loads.large_pressure, 2.0e5)
    np.testing.assert_allclose(loads.small_gas_force, 400.0)
    np.testing.assert_allclose(loads.large_gas_force, 800.0)
    np.testing.assert_allclose(loads.generalized_gas_torque, 0.0)
    assert loads.generalized_net_pressure_torque is None
    assert loads.small_net_pressure_force is None
    assert loads.large_net_pressure_force is None


def test_external_pressure_produces_explicit_net_piston_force(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, _, cycle = create_static_cycle(ideal_gas)
    loads = extract_thermodynamic_loads(
        cycle,
        model,
        PistonFaceAreas(small=2.0e-3, large=4.0e-3),
        PistonExternalPressures(small=1.0e5, large=1.0e5),
    )

    np.testing.assert_allclose(loads.small_net_pressure_force, 200.0)
    np.testing.assert_allclose(loads.large_net_pressure_force, 400.0)
    np.testing.assert_allclose(loads.generalized_net_pressure_torque, 0.0)


def test_piston_area_must_be_explicit_and_positive() -> None:
    with pytest.raises(ValueError, match="positive"):
        PistonFaceAreas(small=0.0, large=1.0e-3)
