"""Simply supported steel IPE 270 beam, 6 m span: 8 kN/m over the whole span plus 20 kN at 2 m.

Run from the repository root:  python examples/beam_example.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mechtools import Beam

E = 210e9          # Pa, structural steel
I = 5790e-8        # m^4, IPE 270 (strong axis)
W = 429e-6         # m^3, elastic section modulus
fy = 235e6         # Pa, S235 yield strength
L = 6.0            # m

beam = Beam(L, E, I).add_uniform_load(8e3).add_point_load(20e3, 2.0)
r = beam.solve()

print(f"R_A = {r.reactions['R_A']/1e3:.2f} kN, R_B = {r.reactions['R_B']/1e3:.2f} kN")
print(f"M_max = {r.max_moment/1e3:.2f} kN·m")
print(f"δ_max = {r.max_deflection*1e3:.2f} mm at x = {r.max_deflection_at:.2f} m "
      f"(L/{L/r.max_deflection:.0f}, limit L/250 -> {'OK' if r.max_deflection <= L/250 else 'FAIL'})")
sigma = abs(r.max_moment) / W
print(f"σ_max = {sigma/1e6:.0f} MPa, utilisation {sigma/fy:.0%} of S235 yield")

ink, accent, fill = "#1b2a33", "#1f6f8b", "#cfe3ea"
fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True, constrained_layout=True)
panels = [
    (r.shear / 1e3, "Shear force V [kN]"),
    (r.moment / 1e3, "Bending moment M [kN·m]"),
    (r.deflection * 1e3, "Deflection δ [mm] (downward +)"),
]
for ax, (y, label) in zip(axes, panels):
    ax.fill_between(r.x, y, color=fill)
    ax.plot(r.x, y, color=accent, lw=2)
    ax.axhline(0, color=ink, lw=0.8)
    ax.set_ylabel(label, color=ink)
    ax.grid(alpha=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[2].invert_yaxis()
axes[1].annotate(f"M_max = {r.max_moment/1e3:.1f} kN·m",
                 (r.x[abs(r.moment).argmax()], r.max_moment / 1e3),
                 xytext=(40, 8), textcoords="offset points", color=ink)
axes[2].annotate(f"δ_max = {r.max_deflection*1e3:.1f} mm",
                 (r.max_deflection_at, r.max_deflection * 1e3),
                 xytext=(-60, 30), textcoords="offset points", color=ink)
axes[2].set_xlabel("x [m]")
fig.suptitle("IPE 270, L = 6 m: 8 kN/m + 20 kN at x = 2 m", color=ink, fontsize=13)
fig.savefig("docs/beam_example.png", dpi=150)
print("saved docs/beam_example.png")
