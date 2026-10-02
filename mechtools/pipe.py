"""Steady incompressible flow in a straight, full, round pipe.

* Reynolds number  Re = ρ V D / μ
* Friction factor  laminar f = 64/Re, turbulent from the Colebrook–White equation
* Pressure drop    Δp = f · (L/D) · ρ V² / 2   (Darcy–Weisbach)

SI units throughout.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

LAMINAR_LIMIT = 2300.0
TURBULENT_LIMIT = 4000.0


def reynolds(velocity: float, diameter: float, density: float, viscosity: float) -> float:
    """Reynolds number from mean velocity [m/s], diameter [m], ρ [kg/m³], μ [Pa·s]."""
    return density * velocity * diameter / viscosity


def haaland(Re: float, rel_roughness: float) -> float:
    """Explicit Haaland (1983) approximation, within ~2 % of Colebrook."""
    return (-1.8 * math.log10((rel_roughness / 3.7) ** 1.11 + 6.9 / Re)) ** -2


def colebrook(Re: float, rel_roughness: float, tol: float = 1e-12, max_iter: int = 100) -> float:
    """Solve 1/√f = -2 log10(ε/3.7D + 2.51/(Re √f)) by fixed-point iteration on 1/√f."""
    x = 1.0 / math.sqrt(haaland(Re, rel_roughness))   # good starting point
    for _ in range(max_iter):
        x_new = -2.0 * math.log10(rel_roughness / 3.7 + 2.51 * x / Re)
        if abs(x_new - x) < tol:
            return 1.0 / x_new**2
        x = x_new
    raise RuntimeError("Colebrook iteration did not converge")


def flow_regime(Re: float) -> str:
    if Re < LAMINAR_LIMIT:
        return "laminar"
    if Re < TURBULENT_LIMIT:
        return "transitional"
    return "turbulent"


def friction_factor(Re: float, rel_roughness: float = 0.0) -> float:
    """Darcy friction factor. Transitional flow uses the turbulent value (conservative)."""
    if Re <= 0:
        raise ValueError("Reynolds number must be positive")
    if Re < LAMINAR_LIMIT:
        return 64.0 / Re
    return colebrook(Re, rel_roughness)


@dataclass
class PipeFlowResult:
    reynolds: float
    regime: str
    friction_factor: float
    velocity: float          # m/s
    pressure_drop: float     # Pa
    head_loss: float         # m of fluid

    def summary(self) -> str:
        return (f"Re = {self.reynolds:,.0f} ({self.regime}), f = {self.friction_factor:.4f}, "
                f"V = {self.velocity:.2f} m/s, Δp = {self.pressure_drop/1000:.2f} kPa, "
                f"h_f = {self.head_loss:.2f} m")


def pressure_drop(flow_rate: float, diameter: float, length: float, density: float = 998.0,
                  viscosity: float = 1.002e-3, roughness: float = 0.0, g: float = 9.81) -> PipeFlowResult:
    """Pressure drop for a volumetric flow rate [m³/s]. Defaults are water at 20 °C."""
    area = math.pi * diameter**2 / 4.0
    V = flow_rate / area
    Re = reynolds(V, diameter, density, viscosity)
    f = friction_factor(Re, roughness / diameter)
    dp = f * (length / diameter) * density * V**2 / 2.0
    return PipeFlowResult(Re, flow_regime(Re), f, V, dp, dp / (density * g))


# Absolute roughness ε of common pipe materials [m] (typical textbook values)
ROUGHNESS = {
    "drawn tubing (copper, glass)": 1.5e-6,
    "PVC / plastic": 1.5e-6,
    "commercial steel": 4.5e-5,
    "galvanised iron": 1.5e-4,
    "cast iron": 2.6e-4,
    "concrete": 1.0e-3,
}
