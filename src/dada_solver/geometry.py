"""Machine gas volumes, independent of the mechanical linkage model."""

from __future__ import annotations

from dataclasses import dataclass
import math


def _require_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive.")


@dataclass(frozen=True, slots=True)
class CylinderVolumeLimits:
    """Minimum and maximum thermodynamic volume of one cylinder in m^3."""

    minimum: float
    maximum: float

    def __post_init__(self) -> None:
        _require_positive("minimum volume", self.minimum)
        _require_positive("maximum volume", self.maximum)
        if self.maximum <= self.minimum:
            raise ValueError("Maximum cylinder volume must exceed minimum volume.")

    @property
    def swept(self) -> float:
        """Return swept volume in m^3."""

        return self.maximum - self.minimum

    def closure_fraction(self, volume: float) -> float:
        """Return Lambda without clipping values outside the nominal limits."""

        return (self.maximum - volume) / self.swept


@dataclass(frozen=True, slots=True)
class MachineVolumes:
    """Cylinder limits and fixed heat-exchanger gas volumes in m^3."""

    small_cylinder: CylinderVolumeLimits
    large_cylinder: CylinderVolumeLimits
    cold_heat_exchanger: float
    hot_heat_exchanger: float

    def __post_init__(self) -> None:
        _require_positive("cold heat-exchanger volume", self.cold_heat_exchanger)
        _require_positive("hot heat-exchanger volume", self.hot_heat_exchanger)


@dataclass(frozen=True, slots=True)
class InstantaneousVolumes:
    """Ordered control-volume values at one instant, all in m^3."""

    small_cylinder: float
    large_cylinder: float
    cold_heat_exchanger: float
    hot_heat_exchanger: float

    def __post_init__(self) -> None:
        for name, value in zip(self.names(), self.as_tuple(), strict=True):
            _require_positive(f"{name} volume", value)

    @staticmethod
    def names() -> tuple[str, ...]:
        return ("S", "L", "C", "H")

    def as_tuple(self) -> tuple[float, float, float, float]:
        return (
            self.small_cylinder,
            self.large_cylinder,
            self.cold_heat_exchanger,
            self.hot_heat_exchanger,
        )

    @property
    def total(self) -> float:
        return sum(self.as_tuple())

