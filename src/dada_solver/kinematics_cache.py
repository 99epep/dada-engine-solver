"""Candidate-local exact reuse for optional combined kinematics providers."""


class ExactAngleKinematics:
    """Reuse only the last combined result; delegate the generic interface.

    The prepared provider must be deterministic for its fixed geometry. This
    adapter never rounds an angle or stores thermodynamic state. Separate scalar
    methods retain their original implementation, including custom providers.
    Create a fresh adapter for each solve; do not share it between candidates.
    """
    def __init__(self, forward):
        self.forward = forward
        self._provider = forward.cylinder_volumes_and_derivatives
        self._angle = None
        self._values = None
        self.calls = self.hits = 0

    def __getattr__(self, name):
        return getattr(self.forward, name)

    def cylinder_volumes_and_derivatives(self, angle):
        self.calls += 1
        if self._values is not None and angle == self._angle:
            self.hits += 1
            return self._values
        # Publish only after successful evaluation, and own an immutable result.
        values = tuple(self._provider(angle))
        self._angle, self._values = angle, values
        return values

    # Explicit protocol members also preserve runtime-checkable Protocol support
    # on Python versions that inspect attributes without invoking __getattr__.
    @property
    def small_volume_limits(self): return self.forward.small_volume_limits

    @property
    def large_volume_limits(self): return self.forward.large_volume_limits

    @property
    def small_physical_stroke(self): return self.forward.small_physical_stroke

    @property
    def large_physical_stroke(self): return self.forward.large_physical_stroke

    def small_cylinder_volume(self, angle): return self.forward.small_cylinder_volume(angle)
    def large_cylinder_volume(self, angle): return self.forward.large_cylinder_volume(angle)
    def small_cylinder_volume_derivative(self, angle): return self.forward.small_cylinder_volume_derivative(angle)
    def large_cylinder_volume_derivative(self, angle): return self.forward.large_cylinder_volume_derivative(angle)
    def breakpoint_angles(self): return getattr(self.forward, 'breakpoint_angles', lambda: ())()

    def snapshot(self):
        return dict(calls=self.calls, hits=self.hits, misses=self.calls-self.hits,
                    hit_fraction=self.hits/self.calls if self.calls else 0.)


def prepare_exact_kinematics(provider):
    """Preserve providers without the optional combined method verbatim."""
    if isinstance(provider, ExactAngleKinematics):
        provider = provider.forward
    if callable(getattr(provider, 'cylinder_volumes_and_derivatives', None)):
        return ExactAngleKinematics(provider)
    return provider
