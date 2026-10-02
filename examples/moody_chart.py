"""Moody chart drawn from the Colebrook–White equation, plus a worked pipe example.

Run from the repository root:  python examples/moody_chart.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mechtools import friction_factor, pressure_drop
from mechtools.pipe import ROUGHNESS

# Worked example: 10 L/s of water, 100 m of DN100 commercial steel
res = pressure_drop(0.010, 0.100, 100.0, roughness=ROUGHNESS["commercial steel"])
print(res.summary())

ink, accent = "#1b2a33", "#1f6f8b"
fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
Re_lam = np.logspace(np.log10(600), np.log10(2300), 50)
ax.plot(Re_lam, 64 / Re_lam, color=ink, lw=1.6)
ax.text(640, 0.034, "laminar\nf = 64/Re", color=ink, fontsize=9)

Re_t = np.logspace(np.log10(4000), 8, 300)
for rr in (0, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2):
    f = [friction_factor(Re, rr) for Re in Re_t]
    ax.plot(Re_t, f, color=accent, lw=1.1, alpha=0.9)
    ax.text(1.15e8, f[-1], "smooth" if rr == 0 else f"{rr:g}", va="center", fontsize=8, color=ink)

ax.axvspan(2300, 4000, color=accent, alpha=0.08)
ax.text(2450, 0.0068, "transition", fontsize=8, color=ink, rotation=90)
ax.plot(res.reynolds, res.friction_factor, "o", color="#c0392b", ms=7, zorder=5)
ax.annotate(f"worked example\nRe = {res.reynolds:,.0f}, f = {res.friction_factor:.4f}",
            (res.reynolds, res.friction_factor), xytext=(15, 25), textcoords="offset points",
            color=ink, fontsize=9, arrowprops=dict(arrowstyle="-", color=ink, lw=0.8))

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(600, 1e8); ax.set_ylim(0.006, 0.1)
ax.set_xlabel("Reynolds number Re"); ax.set_ylabel("Darcy friction factor f")
ax.set_title("Moody chart (Colebrook–White)   ·   labels: relative roughness ε/D", color=ink)
ax.grid(which="both", alpha=0.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.savefig("docs/moody_chart.png", dpi=150)
print("saved docs/moody_chart.png")
