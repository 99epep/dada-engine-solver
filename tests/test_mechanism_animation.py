import math

import pytest

from dada_solver.mechanism_animation import piston_bore


def test_piston_bore_reconstructs_swept_volume() -> None:
    swept_volume = 0.0066
    stroke = 0.20

    bore = piston_bore(swept_volume, stroke)

    assert math.pi * bore**2 * stroke / 4.0 == pytest.approx(swept_volume)


def test_piston_bore_rejects_non_positive_geometry() -> None:
    with pytest.raises(ValueError, match="positive"):
        piston_bore(0.0, 0.2)
