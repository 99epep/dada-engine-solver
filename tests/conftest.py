import pytest

from dada_solver.fluids import CaloricallyPerfectGas


@pytest.fixture
def ideal_gas() -> CaloricallyPerfectGas:
    return CaloricallyPerfectGas(
        gas_constant=287.0,
        heat_capacity_cp=1004.5,
        heat_capacity_cv=717.5,
    )

