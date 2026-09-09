"""Independent periodic cubic volume laws without prescribed cycle phases."""
from dataclasses import dataclass, field, replace
import math
import numpy as np
from scipy.interpolate import CubicSpline
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.kinematics import KinematicConstraintViolation

TAU = 2 * math.pi


@dataclass(frozen=True, slots=True)
class MotionDiagnostics:
    """Analytic cycle extrema and optional signed feasibility margins, in SI."""
    maximum_absolute_first_derivative: float
    maximum_absolute_second_derivative: float
    first_derivative_margin: float | None
    second_derivative_margin: float | None

    @property
    def feasible(self) -> bool:
        return all(x is None or x >= 0 for x in
                   (self.first_derivative_margin, self.second_derivative_margin))


@dataclass(frozen=True, slots=True)
class FreeMotionDefinition:
    """N equally spaced shape controls; physical offset/amplitude live in limits.

    Controls are canonicalized to zero mean and unit norm. For unconstrained
    optimization use N-2 stereographic shape coordinates, rather than N raw
    controls with redundant affine degrees of freedom. No phase is prescribed.
    """
    control_values: tuple[float, ...]
    minimum_volume: float
    maximum_volume: float
    maximum_absolute_first_derivative: float | None = None
    maximum_absolute_second_derivative: float | None = None

    def __post_init__(self):
        CylinderVolumeLimits(self.minimum_volume, self.maximum_volume)
        values = np.asarray(self.control_values, dtype=float)
        if values.ndim != 1 or len(values) < 4 or not np.all(np.isfinite(values)):
            raise ValueError('At least four finite scalar control values are required.')
        # Scale before centering to avoid overflow for large but finite inputs.
        scale = float(np.max(np.abs(values)))
        centered = values / scale if scale else values.copy()
        centered = centered - np.mean(centered)
        norm = float(np.linalg.norm(centered))
        if norm <= 32 * np.finfo(float).eps:
            raise ValueError('Free motion must have a nonzero, numerically resolved shape.')
        canonical = centered / norm
        # Preserve already canonical serialized values bit-for-bit on reconstruction.
        if scale <= 1.0 and abs(float(np.mean(values))) <= 8*np.finfo(float).eps and abs(float(np.linalg.norm(values))-1) <= 8*np.finfo(float).eps:
            canonical = values
        object.__setattr__(self, 'control_values', tuple(float(x) for x in canonical))
        for value in (self.maximum_absolute_first_derivative, self.maximum_absolute_second_derivative):
            if value is not None and (not math.isfinite(value) or value < 0):
                raise ValueError('Derivative limits must be finite and nonnegative.')

    @property
    def limits(self):
        return CylinderVolumeLimits(self.minimum_volume, self.maximum_volume)

    @classmethod
    def from_shape_coordinates(cls, coordinates, minimum_volume, maximum_volume, **limits):
        """Build N controls from N-2 nonredundant coordinates in one sphere chart.

        A finite chart excludes one pole; another chart/control input can cover
        it. This is a parameterization, not a constraint on preferred phases.
        """
        z = np.asarray(coordinates, dtype=float)
        if z.ndim != 1 or len(z) < 2 or not np.all(np.isfinite(z)):
            raise ValueError('At least two finite shape coordinates are required.')
        square = float(z @ z)
        if not math.isfinite(square):
            raise ValueError('Shape coordinates are too large for this chart.')
        sphere = np.r_[2*z, square-1] / (1+square)
        n = len(z)+2
        basis = np.zeros((n, n-1))
        for j in range(n-1):
            basis[:j+1,j] = 1 / math.sqrt((j+1)*(j+2))
            basis[j+1,j] = -(j+1) / math.sqrt((j+1)*(j+2))
        return cls(tuple(basis @ sphere), minimum_volume, maximum_volume, **limits)


@dataclass(frozen=True, slots=True)
class FreeKinematicsConfiguration:
    small: FreeMotionDefinition
    large: FreeMotionDefinition

    def with_volume_limits(self, small, large):
        """Update physical ranges without adding new shape/amplitude variables."""
        return replace(self,
            small=replace(self.small, minimum_volume=small.minimum, maximum_volume=small.maximum),
            large=replace(self.large, minimum_volume=large.minimum, maximum_volume=large.maximum))


@dataclass(frozen=True, slots=True)
class _PeriodicMotion:
    definition: FreeMotionDefinition
    _coefficients: tuple[tuple[float, ...], ...] = field(init=False, repr=False)
    _knots: tuple[float, ...] = field(init=False, repr=False)
    _minimum: float = field(init=False)
    _scale: float = field(init=False)
    diagnostics: MotionDiagnostics = field(init=False)

    def __post_init__(self):
        values = self.definition.control_values
        knots = np.linspace(0, TAU, len(values)+1)
        spline = CubicSpline(knots, (*values, values[0]), bc_type='periodic')
        def candidates(order):
            roots = spline.derivative(order+1).roots(extrapolate=False)
            return np.r_[knots, roots[np.isfinite(roots) & (roots >= 0) & (roots <= TAU)]]
        extrema = spline(candidates(0))
        minimum, maximum = float(np.min(extrema)), float(np.max(extrema))
        scale = self.definition.limits.swept / (maximum-minimum)
        first = float(np.max(np.abs(spline(candidates(1), 1)))) * scale
        second = float(np.max(np.abs(spline(knots, 2)))) * scale
        # Store only immutable polynomial data; evaluate analytic derivatives.
        object.__setattr__(self, '_coefficients', tuple(tuple(float(x) for x in row) for row in spline.c))
        object.__setattr__(self, '_knots', tuple(float(x) for x in knots))
        object.__setattr__(self, '_minimum', minimum)
        object.__setattr__(self, '_scale', scale)
        d = self.definition
        object.__setattr__(self, 'diagnostics', MotionDiagnostics(first, second,
            None if d.maximum_absolute_first_derivative is None else d.maximum_absolute_first_derivative-first,
            None if d.maximum_absolute_second_derivative is None else d.maximum_absolute_second_derivative-second))

    def evaluate(self, theta, order=0):
        angle = np.asarray(theta, dtype=float)
        if not np.all(np.isfinite(angle)):
            raise ValueError('Cycle angle must be finite.')
        wrapped = np.remainder(angle, TAU)
        index = np.searchsorted(self._knots, wrapped, side='right') - 1
        dx = wrapped - np.asarray(self._knots)[index]
        a, b, c, d = np.asarray(self._coefficients)[:, index]
        if order == 0:
            value = ((a*dx+b)*dx+c)*dx+d
        elif order == 1:
            value = (3*a*dx+2*b)*dx+c
        elif order == 2:
            value = 6*a*dx+2*b
        else:
            raise ValueError('Only volume and its first two derivatives are supported.')
        result = ((value-self._minimum)*self._scale+self.definition.minimum_volume
                  if order == 0 else value*self._scale)
        return float(result) if result.ndim == 0 else result


@dataclass(frozen=True, slots=True)
class FreeKinematics:
    """C2 periodic motion on theta in [0, 2*pi), extended by periodic wrapping.

    Definitions use the study-angle convention, as existing backend families
    do. The factory reverses it once for negative-speed motor operation.
    """
    configuration: FreeKinematicsConfiguration
    _small: _PeriodicMotion = field(init=False, repr=False)
    _large: _PeriodicMotion = field(init=False, repr=False)

    def __post_init__(self):
        object.__setattr__(self, '_small', _PeriodicMotion(self.configuration.small))
        object.__setattr__(self, '_large', _PeriodicMotion(self.configuration.large))

    @property
    def small_volume_limits(self):
        return self.configuration.small.limits

    @property
    def large_volume_limits(self):
        return self.configuration.large.limits

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    @property
    def diagnostics(self):
        return (self._small.diagnostics, self._large.diagnostics)

    def require_feasible(self):
        if not all(d.feasible for d in self.diagnostics):
            raise KinematicConstraintViolation('Free kinematics violates configured derivative limits; inspect diagnostics.', self.diagnostics)

    def small_cylinder_volume(self, theta):
        return self._small.evaluate(theta)

    def large_cylinder_volume(self, theta):
        return self._large.evaluate(theta)

    def small_cylinder_volume_derivative(self, theta):
        return self._small.evaluate(theta, 1)

    def large_cylinder_volume_derivative(self, theta):
        return self._large.evaluate(theta, 1)

    def small_cylinder_volume_second_derivative(self, theta):
        return self._small.evaluate(theta, 2)

    def large_cylinder_volume_second_derivative(self, theta):
        return self._large.evaluate(theta, 2)

    def cylinder_volumes_and_derivatives(self, theta):
        return (self.small_cylinder_volume(theta), self.large_cylinder_volume(theta),
                self.small_cylinder_volume_derivative(theta), self.large_cylinder_volume_derivative(theta))

    def breakpoint_angles(self):
        return ()
