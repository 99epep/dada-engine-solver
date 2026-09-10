"""Incremental global exploration; no local solver or objective penalties."""
from typing import Protocol
import warnings
from scipy.stats import qmc


class SearchStrategy(Protocol):
    @property
    def index(self) -> int: ...
    def next_point(self) -> tuple[float, ...]: ...
    def state(self) -> dict: ...


class SobolStrategy:
    def __init__(self, dimension, *, seed=0, scramble=True, index=0):
        if dimension < 1 or index < 0:
            raise ValueError('Dimension must be positive and sequence index nonnegative.')
        self.dimension, self.seed, self.scramble = dimension, seed, scramble
        self._engine = qmc.Sobol(d=dimension, scramble=scramble, seed=seed)
        if index:
            self._engine.fast_forward(index)

    @property
    def index(self):
        return self._engine.num_generated

    def next_point(self):
        # Arbitrary budget endpoints are allowed; balanced Sobol prefixes occur at 2**m.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='The balance properties of Sobol')
            return tuple(float(x) for x in self._engine.random(1)[0])

    def state(self):
        return dict(type='sobol', dimension=self.dimension, seed=self.seed,
                    scramble=self.scramble, index=self.index)
