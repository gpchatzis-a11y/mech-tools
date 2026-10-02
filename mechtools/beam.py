"""Euler–Bernoulli beam analysis: reactions, shear, bending moment and deflection.

Supports
--------
* ``"simply_supported"`` – pin at x = 0, roller at x = L
* ``"cantilever"``       – fixed at x = 0, free at x = L

Sign convention
---------------
Downward loads are positive numbers. Sagging bending moment is positive.
Deflection is positive *downward*, so a loaded beam shows positive values.
Units are whatever you feed in, as long as they are consistent (SI: N, m, Pa).
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass(frozen=True)
class PointLoad:
    P: float   # magnitude (downward positive) [N]
    a: float   # position from the left end [m]


@dataclass(frozen=True)
class UniformLoad:
    w: float          # intensity (downward positive) [N/m]
    start: float = 0.0
    end: float | None = None   # None = to the end of the beam


@dataclass
class BeamResult:
    x: np.ndarray
    shear: np.ndarray
    moment: np.ndarray
    deflection: np.ndarray
    reactions: dict

    @property
    def max_deflection(self) -> float:
        return float(np.max(np.abs(self.deflection)))

    @property
    def max_deflection_at(self) -> float:
        return float(self.x[np.argmax(np.abs(self.deflection))])

    @property
    def max_moment(self) -> float:
        return float(self.moment[np.argmax(np.abs(self.moment))])


@dataclass
class Beam:
    length: float                    # L [m]
    E: float                         # Young's modulus [Pa]
    I: float                         # second moment of area [m^4]
    support: str = "simply_supported"
    loads: list = field(default_factory=list)

    def __post_init__(self):
        if self.support not in ("simply_supported", "cantilever"):
            raise ValueError("support must be 'simply_supported' or 'cantilever'")
        if self.length <= 0 or self.E <= 0 or self.I <= 0:
            raise ValueError("length, E and I must be positive")

    # -- building the load case ---------------------------------------------
    def add_point_load(self, P: float, a: float) -> "Beam":
        if not 0 <= a <= self.length:
            raise ValueError("point load position must lie on the beam")
        self.loads.append(PointLoad(P, a))
        return self

    def add_uniform_load(self, w: float, start: float = 0.0, end: float | None = None) -> "Beam":
        end = self.length if end is None else end
        if not 0 <= start < end <= self.length:
            raise ValueError("uniform load must satisfy 0 <= start < end <= L")
        self.loads.append(UniformLoad(w, start, end))
        return self

    # -- statics -------------------------------------------------------------
    def _resultants(self):
        """Each load as (force, centroid) so equilibrium is one sum."""
        out = []
        for ld in self.loads:
            if isinstance(ld, PointLoad):
                out.append((ld.P, ld.a))
            else:
                end = self.length if ld.end is None else ld.end
                out.append((ld.w * (end - ld.start), 0.5 * (ld.start + end)))
        return out

    def reactions(self) -> dict:
        res = self._resultants()
        total = sum(F for F, _ in res)
        if self.support == "simply_supported":
            RB = sum(F * x for F, x in res) / self.length
            return {"R_A": total - RB, "R_B": RB}
        # cantilever, fixed at x = 0: vertical reaction and (hogging) fixed-end moment
        return {"R_A": total, "M_A": -sum(F * x for F, x in res)}

    # -- internal forces at a section ---------------------------------------
    def _shear_moment(self, x: np.ndarray):
        r = self.reactions()
        V = np.full_like(x, r["R_A"], dtype=float)
        M = r["R_A"] * x + r.get("M_A", 0.0)
        for ld in self.loads:
            if isinstance(ld, PointLoad):
                left = x > ld.a
                V -= np.where(left, ld.P, 0.0)
                M -= np.where(left, ld.P * (x - ld.a), 0.0)
            else:
                end = self.length if ld.end is None else ld.end
                covered = np.clip(x - ld.start, 0.0, end - ld.start)   # loaded length left of x
                F = ld.w * covered
                arm = x - (ld.start + covered / 2)
                V -= F
                M -= F * arm
        return V, M

    # -- full solution -------------------------------------------------------
    def solve(self, n: int = 2001) -> BeamResult:
        x = np.linspace(0.0, self.length, n)
        V, M = self._shear_moment(x)
        kappa = M / (self.E * self.I)   # curvature, EI·v'' = -M with v downward

        # integrate twice with the trapezoidal rule
        dx = np.diff(x)
        theta = np.concatenate(([0.0], np.cumsum(0.5 * (kappa[1:] + kappa[:-1]) * dx)))
        v = np.concatenate(([0.0], np.cumsum(0.5 * (theta[1:] + theta[:-1]) * dx)))
        v = -v  # downward positive

        if self.support == "simply_supported":
            # add a linear term so that v(0) = v(L) = 0
            v = v - v[-1] * x / self.length
        # cantilever: v(0) = 0 and v'(0) = 0 are already satisfied

        return BeamResult(x=x, shear=V, moment=M, deflection=v, reactions=self.reactions())


def rectangular_I(b: float, h: float) -> float:
    """Second moment of area of a solid b × h rectangle about its neutral axis."""
    return b * h**3 / 12.0


def circular_I(d: float) -> float:
    """Second moment of area of a solid circular section of diameter d."""
    return np.pi * d**4 / 64.0
