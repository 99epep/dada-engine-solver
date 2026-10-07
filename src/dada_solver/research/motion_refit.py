"""Scientific refit request boundary; no approximate initializer is supplied."""
from dataclasses import dataclass
from .motion_target import MotionTarget


@dataclass(frozen=True)
class MotionRefitRequest:
    target: MotionTarget
    destination_family: str = 'free_spline'
    points_per_piston: int = 15
    position_role: str = 'primary_synthesis_reference'
    acceleration_role: str = 'diagnostic_only'

    def __post_init__(self):
        if not isinstance(self.target, MotionTarget):
            raise ValueError('Refit requires a MotionTarget.')
        if self.destination_family != 'free_spline' or type(self.points_per_piston) is not int or self.points_per_piston != 15:
            raise ValueError('The prepared refit contract is free_spline with 15 points per piston.')
        if (self.position_role, self.acceleration_role) != ('primary_synthesis_reference', 'diagnostic_only'):
            raise ValueError('Refit uses position as reference and acceleration as diagnostic only.')

    def execute(self):
        raise NotImplementedError('Motion refit is not implemented; the feature-aware free_spline initializer with 15 points per piston is a separate scientific step.')
