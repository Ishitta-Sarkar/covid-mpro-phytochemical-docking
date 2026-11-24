"""Coordinate geometry used for the search box and interaction criteria."""

from __future__ import annotations

import math
from typing import Iterable, Sequence

from mpro_docking.constants import (
    BOX_PADDING,
    CENTER_DECIMALS,
    MIN_BOX_SIZE,
)

Vec = Sequence[float]


def distance(a: Vec, b: Vec) -> float:
    return math.sqrt(sum((float(a[i]) - float(b[i])) ** 2 for i in range(3)))


def angle_degrees(a: Vec, b: Vec, c: Vec) -> float:
    """Return the angle at ``b`` formed by points a-b-c, in degrees."""

    v1 = [float(a[i]) - float(b[i]) for i in range(3)]
    v2 = [float(c[i]) - float(b[i]) for i in range(3)]
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0.0 or n2 == 0.0:
        raise ValueError("angle is undefined for a zero-length vector")
    cos = sum(v1[i] * v2[i] for i in range(3)) / (n1 * n2)
    cos = max(-1.0, min(1.0, cos))
    return math.degrees(math.acos(cos))


def centroid(coords: Iterable[Vec]) -> tuple[float, float, float]:
    pts = [tuple(map(float, c)) for c in coords]
    if not pts:
        raise ValueError("centroid of an empty coordinate list")
    n = float(len(pts))
    return tuple(sum(p[i] for p in pts) / n for i in range(3))  # type: ignore[return-value]


def cubic_search_box(
    coords: Iterable[Vec],
    padding: float = BOX_PADDING,
    minimum: float = MIN_BOX_SIZE,
    decimals: int = CENTER_DECIMALS,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Return a cubic Vina box that contains ``coords``.

    The centre is the heavy-atom centroid, rounded to ``decimals`` places.
    The edge is the smallest even number of angstroms that is at least
    ``minimum`` and at least ``2 * (max absolute deviation + padding)``.
    An even edge keeps the half-width on the 1 Å grid Vina users expect.
    """

    pts = [tuple(map(float, c)) for c in coords]
    center = centroid(pts)
    half = max(max(abs(p[i] - center[i]) for i in range(3)) for p in pts)
    needed = 2.0 * (half + padding)
    size = max(float(minimum), 2.0 * math.ceil(needed / 2.0))
    rounded = tuple(round(v, decimals) for v in center)
    return rounded, (size, size, size)  # type: ignore[return-value]


def is_hydrogen_bond(
    donor: Vec,
    hydrogen: Vec,
    acceptor: Vec,
    max_ha: float = 2.5,
    min_angle: float = 120.0,
) -> tuple[bool, float, float]:
    """Baker-Hubbard test. Returns (accepted, H...A distance, D-H...A angle)."""

    ha = distance(hydrogen, acceptor)
    ang = angle_degrees(donor, hydrogen, acceptor)
    return ha <= max_ha and ang >= min_angle, ha, ang
