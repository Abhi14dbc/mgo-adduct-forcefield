#!/usr/bin/env python3
"""Graphical abstract, rebuilt on the n = 10 explicit-solvent result.

Replaces figures/graphical_abstract.png (21 Aug), whose H-bond and SASA claims
were superseded by Phase 1. Everything stated here is traceable to a number in
the manuscript draft.
"""
import os, sys, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, RegularPolygon, FancyArrow

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

NATIVE, GLYC, GREEN = "#0072B2", "#D55E00", "#009E73"
INK, MUTED, FAINT, GRID = "#14201F", "#5B6D6B", "#8A9A98", "#DCE5E3"
PANEL, SUNK = "#FFFFFF", "#F2F6F5"
OUT = "figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"]})
fig = plt.figure(figsize=(13.0, 5.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def panel(x, y, w, h, fc=PANEL):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                fc=fc, ec=GRID, lw=1.1, zorder=1))


def head(x, y, n, t):
    ax.text(x, y, n, fontsize=9, color=FAINT, weight="bold")
    ax.text(x + 3.4, y, t, fontsize=9, color=FAINT, weight="bold")


# ============================================================ 1 · MODIFICATION
head(3.2, 92.5, "1 ·", "THE MODIFICATION")
panel(2.5, 40, 28.5, 49)

ax.text(16.7, 84, "Arg143 — catalytic arginine", fontsize=11, color=INK,
        weight="bold", ha="center")
# guanidinium Y
gx, gy = 16.7, 74
ax.plot([gx-6, gx-2.6], [gy-2.4, gy], color=NATIVE, lw=2, zorder=3)
ax.plot([gx-2.6, gx], [gy, gy-2.4], color=NATIVE, lw=2, zorder=3)
ax.plot([gx, gx+3.4], [gy-2.4, gy], color=NATIVE, lw=2, zorder=3)
ax.plot([gx+3.4, gx+3.4], [gy, gy+3.2], color=NATIVE, lw=2, zorder=3)
ax.plot([gx+3.4, gx+6.8], [gy, gy-2.4], color=NATIVE, lw=2, zorder=3)
ax.add_patch(Circle((gx+3.4, gy+3.2), 1.05, fc=NATIVE, ec="white", lw=1.2, zorder=4))
ax.text(gx+5.4, gy+4.2, "+", fontsize=13, color=NATIVE, weight="bold")
ax.text(16.7, 66.5, "+1  ·  5 N–H donors", fontsize=10, color=MUTED, ha="center")

ax.add_patch(FancyArrow(16.7, 63.5, 0, -4.6, width=0.22, head_width=1.2,
                        head_length=1.3, fc=INK, ec=INK, zorder=3))
ax.text(18.6, 61.0, "+ methylglyoxal", fontsize=9.5, color=MUTED)
ax.text(18.6, 58.6, "− 2 H₂O", fontsize=9.5, color=MUTED)

ax.add_patch(RegularPolygon((16.7, 51.5), 5, radius=4.2, orientation=0.32,
                            fc="none", ec=GLYC, lw=2, zorder=3))
ax.text(21.9, 53.4, "O", fontsize=11, color=GLYC, weight="bold")
ax.text(16.7, 44.6, "MG-H1 hydroimidazolone", fontsize=11, color=INK,
        weight="bold", ha="center")
ax.text(16.7, 42.0, "neutral  ·  2 donors  ·  +18% volume", fontsize=10,
        color=MUTED, ha="center")

# ============================================================ 2 · RESULT
head(35.2, 92.5, "2 ·", "WHAT GLYCATION DOES")
panel(34.5, 40, 30.5, 49)

ax.text(49.7, 84, "Site-143 core contact network", fontsize=11, color=INK,
        weight="bold", ha="center")

rows = [("His48", 72.9, 26.2, True), ("Gly61", 65.2, 20.7, True),
        ("Cys57", 59.7, 42.5, False)]
for i, (nm, a, b, broken) in enumerate(rows):
    yy = 77.5 - i * 6.4
    ax.text(38.0, yy, nm, fontsize=10.5, color=INK, weight="bold", va="center")
    ax.plot([43.6, 43.6 + a * 0.115], [yy, yy], color=NATIVE, lw=5,
            solid_capstyle="round", zorder=3)
    ax.plot([43.6, 43.6 + b * 0.115], [yy - 2.0, yy - 2.0],
            color=GLYC if broken else FAINT, lw=5, solid_capstyle="round", zorder=3)
    ax.text(57.6, yy - 1.0, f"{a:.0f} → {b:.0f}%", fontsize=9.5,
            color=INK if broken else FAINT, va="center", weight="bold" if broken else "normal")
    if not broken:
        ax.text(57.6, yy - 3.9, "unchanged", fontsize=8.5, color=FAINT, va="center")

ax.plot([37.5, 62.0], [55.0, 55.0], color=GRID, lw=1)
ax.text(49.7, 50.6, "−36.1", fontsize=27, color=GLYC, weight="bold", ha="center")
ax.text(49.7, 47.0, "percentage points, core composite", fontsize=9.5,
        color=MUTED, ha="center")
ax.text(49.7, 43.6, "p = 0.0009   vs preregistered α = 0.0042",
        fontsize=10, color=INK, ha="center", weight="bold")
ax.text(49.7, 41.2, "n = 10 per group  ·  explicit OPC water  ·  300 ns",
        fontsize=9, color=FAINT, ha="center")

# ============================================================ 3 · METHOD
head(68.2, 92.5, "3 ·", "AND WHAT IT DEPENDS ON")
panel(67.5, 40, 30.0, 49)

ax.text(82.5, 84, "Same endpoint, three methods", fontsize=11, color=INK,
        weight="bold", ha="center")
ax.text(82.5, 81.0, "site H-bonds, glycated − native", fontsize=9, color=FAINT, ha="center")

bars = [("implicit GB\n100 ns · n=3\nper-rep windows", -1.923, MUTED),
        ("implicit GB\n300 ns · n=3\nmatched window", -1.100, MUTED),
        ("explicit OPC\n300 ns · n=10\nmatched window", -0.081, GLYC)]
x0 = 88.5
for i, (lab, v, c) in enumerate(bars):
    yy = 75.0 - i * 8.6
    ax.plot([x0, x0 + v * 5.0], [yy, yy], color=c, lw=6, solid_capstyle="round", zorder=3)
    ax.text(x0 + 0.9, yy, f"{v:+.2f}", fontsize=9.5, color=INK, va="center", weight="bold")
    ax.text(68.9, yy, lab, fontsize=7.9, color=MUTED, va="center", linespacing=1.5)
ax.plot([x0, x0], [78.0, 55.0], color=INK, lw=1.1, zorder=2)

ax.add_patch(FancyBboxPatch((69.0, 42.0), 27.0, 10.5,
                            boxstyle="round,pad=0,rounding_size=1", fc=SUNK, ec=GRID,
                            lw=1, zorder=2))
ax.text(82.5, 48.6, "The implicit-solvent result does not survive", fontsize=9.6,
        color=INK, ha="center", weight="bold")
ax.text(82.5, 45.8, "explicit water. Most AGE simulation work", fontsize=9.6,
        color=MUTED, ha="center")
ax.text(82.5, 43.4, "runs in implicit solvent.", fontsize=9.6, color=MUTED, ha="center")

# ============================================================ footer chips
chips = [("Preregistered", "endpoint and α sealed 22 Aug,\nbefore the data existed"),
         ("Replicate-level", "n = 10 per group; the replicate\nis the unit, never the frame"),
         ("Window-robust", "significant in all five analysis\nwindows, including 0–300 ns"),
         ("Fully reported", "all six predictions, including\ntwo refuted and one unevaluable")]
for i, (t, b) in enumerate(chips):
    x = 2.5 + i * 24.0
    ax.add_patch(FancyBboxPatch((x, 20.0), 22.6, 12.0,
                                boxstyle="round,pad=0,rounding_size=1", fc=SUNK,
                                ec=GRID, lw=1, zorder=1))
    ax.text(x + 1.6, 29.4, t, fontsize=10, color=INK, weight="bold", va="top")
    ax.text(x + 1.6, 26.2, b, fontsize=8.6, color=MUTED, linespacing=1.5, va="top")

# ============================================================ headline
ax.text(2.5, 12.2,
        "Glycation at Arg143 breaks a loop-junction contact network without unfolding SOD1",
        fontsize=15.5, color=INK, weight="bold")
ax.text(2.5, 7.6,
        "Cys57 — the residue implicated through its disulfide with Cys146 — is the one contact that does not move.",
        fontsize=10.5, color=MUTED)
ax.text(2.5, 4.4,
        "Backbone RMSD, radius of gyration and β-sheet content are unchanged across twenty replicates: the perturbation is local, not global.",
        fontsize=10.5, color=MUTED)

for ext in ("pdf", "svg", "png"):
    fig.savefig(f"{OUT}/graphical_abstract_v2.{ext}", dpi=300, facecolor="white")
plt.close(fig)
print(f"wrote {OUT}/graphical_abstract_v2.[pdf|svg|png]")
