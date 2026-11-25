"""Transparent engineering calculations for hydraulic jumps.

The module deliberately keeps the assumptions visible: prismatic channel,
hydrostatic pressure in the far sections, and momentum conservation across
the jump. It is intended for preliminary design and checking, not as a
replacement for a project standard or physical model.
"""

from dataclasses import dataclass
from math import sqrt
from typing import Literal

Section = Literal["rectangular", "trapezoidal"]


@dataclass(frozen=True)
class JumpResult:
    q: float
    y1: float
    y2: float
    v1: float
    v2: float
    froude1: float
    froude2: float
    critical_depth: float
    energy_loss: float
    specific_energy1: float
    specific_energy2: float
    jump_length: float
    regime: str
    basin_recommendation: str
    warnings: tuple[str, ...]


def _geometry(y: float, width: float, side_slope: float, section: Section):
    if section == "rectangular":
        area = width * y
        top_width = width
        centroid_depth = y / 2.0
    else:
        # side_slope is horizontal:vertical (z); width is bottom width.
        area = y * (width + side_slope * y)
        top_width = width + 2.0 * side_slope * y
        centroid_depth = y * (width / 2.0 + side_slope * y / 3.0) / (width + side_slope * y)
    return area, top_width, centroid_depth


def _momentum(q: float, y: float, width: float, side_slope: float, section: Section) -> float:
    area, _, centroid_depth = _geometry(y, width, side_slope, section)
    return q * q / (9.80665 * area) + area * centroid_depth


def _bisect(fn, lo: float, hi: float, iterations: int = 120) -> float:
    flo, fhi = fn(lo), fn(hi)
    while flo * fhi > 0 and hi < 1e8:
        hi *= 2.0
        fhi = fn(hi)
    if flo * fhi > 0:
        raise ValueError("Не удалось найти сопряжённую глубину в заданном диапазоне.")
    for _ in range(iterations):
        mid = (lo + hi) / 2.0
        fmid = fn(mid)
        if abs(fmid) < 1e-11:
            return mid
        if flo * fmid <= 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return (lo + hi) / 2.0


def calculate_jump(
    q: float,
    y1: float,
    width: float,
    *,
    section: Section = "rectangular",
    side_slope: float = 0.0,
    tailwater: float | None = None,
) -> JumpResult:
    if q <= 0 or y1 <= 0 or width <= 0:
        raise ValueError("Q, y1 и ширина должны быть положительными.")
    if section == "trapezoidal" and side_slope < 0:
        raise ValueError("Уклон откоса z должен быть неотрицательным.")

    g = 9.80665
    a1, t1, _ = _geometry(y1, width, side_slope, section)
    v1 = q / a1
    hydraulic_depth = a1 / t1
    fr1 = v1 / sqrt(g * hydraulic_depth)
    if fr1 <= 1.0:
        raise ValueError("Прыжок возможен только при сверхкритическом входном потоке: Fr1 > 1.")

    m1 = _momentum(q, y1, width, side_slope, section)
    y2 = _bisect(lambda y: _momentum(q, y, width, side_slope, section) - m1, y1 * (1 + 1e-8), max(y1 * 2, 1.0))
    a2, t2, _ = _geometry(y2, width, side_slope, section)
    v2 = q / a2
    fr2 = v2 / sqrt(g * (a2 / t2))
    e1 = y1 + v1 * v1 / (2 * g)
    e2 = y2 + v2 * v2 / (2 * g)
    loss = e1 - e2
    critical = _bisect(
        lambda y: q * q * _geometry(y, width, side_slope, section)[1]
        / (g * (_geometry(y, width, side_slope, section)[0] ** 3))
        - 1.0,
        1e-8,
        max(y2 * 2, 1.0),
    )
    jump_length = y2 * (4.5 if fr1 < 4.5 else 6.0)

    if fr1 < 1.7:
        regime = "переходная область; устойчивость прыжка низкая"
        basin = "Свободный прыжок применять осторожно; проверить физической моделью или 2D-расчётом."
    elif fr1 < 4.5:
        regime = "слабый/переходный прыжок"
        basin = "Ориентировочно USBR Type IV или SAF; проверить диапазон стандарта."
    else:
        regime = "развитый турбулентный прыжок"
        basin = "Ориентировочно USBR Type II/III или SAF; окончательный выбор по стандарту сооружения."

    warnings = []
    if tailwater is not None and tailwater < y2:
        warnings.append("Подпор ниже сопряжённой глубины: прыжок может сместиться вниз по потоку.")
    if tailwater is not None and tailwater > 1.2 * y2:
        warnings.append("Подпор существенно выше сопряжённой глубины: возможен затопленный прыжок.")
    if fr1 < 1.7 or fr1 > 17:
        warnings.append("Число Фруда вне типового диапазона применения распространённых бассейнов SAF/USBR.")
    return JumpResult(q, y1, y2, v1, v2, fr1, fr2, critical, loss, e1, e2, jump_length, regime, basin, tuple(warnings))
