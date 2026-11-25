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


def test_zero_side_slope_matches_rectangular():
    rectangular = calculate_jump(10, 0.45, 3)
    trapezoidal = calculate_jump(10, 0.45, 3, section="trapezoidal", side_slope=0)
    assert math.isclose(rectangular.y2, trapezoidal.y2, rel_tol=1e-9)


def test_trapezoidal_momentum_balance():
    result = calculate_jump(10, 0.45, 3, section="trapezoidal", side_slope=1.5)
    def force(y):
        return 100 / (9.80665 * y * (3 + 1.5 * y)) + 3 * y**2 / 2 + 1.5 * y**3 / 3
    assert math.isclose(force(result.y1), force(result.y2), rel_tol=1e-9)
    assert result.froude2 < 1
