import math

from hydraulic_jump.core import calculate_jump


def test_rectangular_jump_matches_belanger_relation():
    result = calculate_jump(q=10.0, y1=0.45, width=3.0)
    expected_ratio = 0.5 * (math.sqrt(1 + 8 * result.froude1**2) - 1)
    assert math.isclose(result.y2 / result.y1, expected_ratio, rel_tol=1e-8)
    assert result.y2 > result.y1
    assert result.energy_loss > 0


def test_trapezoidal_jump_is_supported():
    result = calculate_jump(q=10.0, y1=0.45, width=3.0, section="trapezoidal", side_slope=1.5)
    assert result.y2 > result.y1
