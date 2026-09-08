#!/usr/bin/env python3
"""Publication figures for the MGO-SOD1 manuscript, from the analysed numbers.

Palette: Okabe-Ito blue (#0072B2) and vermillion (#D55E00) -- the published
colour-universal-design pair. CVD separation is computed below rather than
assumed. Outputs PDF (vector, for submission), SVG (for editing) and PNG.
"""
import io, sys, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

NATIVE, GLYC = "#0072B2", "#D55E00"
INK, MUTED, GRID = "#14201F", "#5B6D6B", "#DCE5E3"
OUT = os.environ.get("FIGDIR", "figures")
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------- CVD check
def _lin(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _oklab(rgb):
    r, g, b = rgb
    l = 0.4122214708*r + 0.5363325363*g + 0.0514459929*b
    m = 0.2119034982*r + 0.6806995451*g + 0.1073969566*b
    s = 0.0883024619*r + 0.2817188376*g + 0.6299787005*b
    l_, m_, s_ = np.cbrt([l, m, s])
    return np.array([0.2104542553*l_ + 0.7936177850*m_ - 0.0040720468*s_,
                     1.9779984951*l_ - 2.4285922050*m_ + 0.4505937099*s_,
                     0.0259040371*l_ + 0.7827717662*m_ - 0.8086757660*s_])


def _cvd(rgb, kind):
    R, G, B = rgb
    L = 17.8824*R + 43.5161*G + 4.11935*B
    M = 3.45565*R + 27.1554*G + 3.86714*B
    S = 0.0299566*R + 0.184309*G + 1.46709*B
    if kind == "deuter":
        M = 0.494207*L + 1.24827*S
    else:
        L = 2.02344*M - 2.52581*S
    inv = np.array([[0.080944, -0.130504, 0.116721],
                    [-0.010248, 0.054019, -0.113615],
                    [-0.000365, -0.004122, 0.693513]])
    return np.clip(inv @ np.array([L, M, S]), 0, 1)


def check(h1, h2):
    rgb = [np.array([_lin(int(h[i:i+2], 16)/255) for i in (1, 3, 5)]) for h in (h1, h2)]
    out = {}
    for kind in ("normal", "deuter", "protan"):
        pair = rgb if kind == "normal" else [_cvd(c, kind) for c in rgb]
        d = np.linalg.norm(_oklab(pair[0]) - _oklab(pair[1])) * 100
        out[kind] = d
    return out


d = check(NATIVE, GLYC)
print("PALETTE CHECK  " + NATIVE + " vs " + GLYC)
for k, v in d.items():
    tgt = 15 if k == "normal" else 8
    print(f"  {k:8s} deltaE = {v:6.2f}   target >= {tgt:2d}   "
          f"{'PASS' if v >= tgt else 'FAIL'}")
print()

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9.5, "axes.titlesize": 10.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{OUT}/{name}.{ext}")
    plt.close(fig)
    print(f"  wrote {OUT}/{name}.[pdf|svg|png]")


# ---------------------------------------------------------------- F1 strip
nat = [9.9, 21.6, 41.9, 45.9, 47.2, 75.5, 92.8, 96.4, 97.3, 99.3]
gly = [0.0, 2.7, 4.3, 5.2, 6.7, 6.8, 8.5, 12.9, 42.4, 100.0]
fig, ax = plt.subplots(figsize=(3.6, 4.2))
rng = np.random.default_rng(7)
for i, (v, c, lab) in enumerate([(nat, NATIVE, "Native"), (gly, GLYC, "Glycated")]):
    x = i + rng.uniform(-0.09, 0.09, len(v))
    ax.scatter(x, v, s=46, facecolor=c, edgecolor="white", linewidth=1.1, zorder=3)
    ax.plot([i-0.26, i+0.26], [np.mean(v)]*2, color=INK, lw=2.2, zorder=4)
    ax.text(i + 0.31, np.mean(v), f"{np.mean(v):.1f}%", ha="left", va="center",
            fontsize=9, color=INK, weight="bold")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Native", "Glycated"])
ax.set_ylabel("Site-143–Gly61 contact occupancy (%)")
ax.set_ylim(-4, 106); ax.set_xlim(-0.5, 1.85)
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
ax.set_axisbelow(True)
ax.set_title("Each point is one 300 ns replicate", fontsize=9, color=MUTED, pad=9)
save(fig, "fig1_replicate_strip")

# ---------------------------------------------------------------- F2 forest
win = ["0–150", "50–200", "100–250", "150–300", "0–300"]
dlt = [-52.39, -48.79, -44.20, -36.14, -44.26]
lo = [-62.76, -61.53, -57.93, -51.33, -55.72]
hi = [-39.79, -34.89, -29.49, -19.83, -31.50]
pv = ["<0.0001", "0.0001", "0.0002", "0.0009", "0.0001"]
fig, ax = plt.subplots(figsize=(5.4, 2.9))
y = np.arange(len(win))[::-1]
for i, (yy, dd, l, h) in enumerate(zip(y, dlt, lo, hi)):
    rep = (win[i] == "150–300")
    ax.plot([l, h], [yy, yy], color=GLYC if rep else MUTED, lw=2.6 if rep else 1.6,
            solid_capstyle="round", zorder=3)
    ax.scatter([dd], [yy], s=70 if rep else 44, facecolor=GLYC if rep else MUTED,
               edgecolor="white", linewidth=1.1, zorder=4)
    ax.text(3, yy, pv[i], va="center", fontsize=8.5, color=INK if rep else MUTED,
            weight="bold" if rep else "normal")
ax.axvline(0, color=INK, lw=1.0, zorder=2)
ax.text(3, len(win)-0.62, "p", fontsize=8.5, color=MUTED, style="italic")
ax.set_yticks(y); ax.set_yticklabels([w + ("  ◀ reported" if w == "150–300" else "") for w in win])
ax.set_xlabel("Δ core contact occupancy, glycated − native (percentage points)")
ax.set_xlim(-70, 14); ax.set_ylim(-0.6, len(win)-0.4)
ax.grid(axis="x", color=GRID, lw=0.7, zorder=0); ax.set_axisbelow(True)
ax.set_title("Significant in every window, including the full trajectory",
             fontsize=9, color=MUTED, pad=9)
save(fig, "fig2_window_forest")

# ---------------------------------------------------------------- F3 decay
blk = [25, 75, 125, 175, 225, 275]
nb = [94.5, 90.9, 95.5, 81.4, 59.4, 47.5]
gb = [28.6, 27.8, 26.4, 21.0, 14.1, 21.7]
fig, ax = plt.subplots(figsize=(5.0, 3.2))
ax.axvspan(150, 300, color=GRID, alpha=0.45, lw=0, zorder=0)
ax.text(225, 103, "analysis window", ha="center", fontsize=8, color=MUTED)
for v, c, lab, m in [(nb, NATIVE, "Native", "o"), (gb, GLYC, "Glycated", "s")]:
    ax.plot(blk, v, color=c, lw=2, marker=m, ms=7, mec="white", mew=1.1,
            label=lab, zorder=3)
ax.set_xlabel("Simulated time (ns, block midpoint)")
ax.set_ylabel("Gly61 contact occupancy (%)")
ax.set_ylim(0, 108); ax.set_xlim(0, 300)
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.grid(color=GRID, lw=0.7, zorder=0); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.62, 0.62))
ax.set_title("The native group is still relaxing at 300 ns", fontsize=9, color=MUTED, pad=9)
save(fig, "fig3_block_decay")

# ---------------------------------------------------------------- F4 H-bond
lbl = ["Implicit GB-OBC2\n100 ns, n=3\nper-replicate windows",
       "Implicit GB-OBC2\n300 ns, n=3\nmatched window",
       "Explicit OPC\n300 ns, n=10\nmatched window"]
dd = [-1.923, -1.100, -0.081]
ll = [-2.438, -1.896, -0.412]
hh = [-1.483, -0.276, 0.303]
fig, ax = plt.subplots(figsize=(5.2, 3.0))
y = np.arange(3)[::-1]
cols = [MUTED, MUTED, GLYC]
for yy, d_, l, h, c in zip(y, dd, ll, hh, cols):
    ax.plot([l, h], [yy, yy], color=c, lw=2.4, solid_capstyle="round", zorder=3)
    ax.scatter([d_], [yy], s=64, facecolor=c, edgecolor="white", linewidth=1.1, zorder=4)
    ax.text(d_, yy + 0.22, f"{d_:+.2f}", ha="center", fontsize=8.5, color=INK)
ax.axvline(0, color=INK, lw=1.0, zorder=2)
ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=8.2)
ax.set_xlabel("Δ site-143 protein H-bonds, glycated − native")
ax.set_xlim(-2.7, 0.7); ax.set_ylim(-0.6, 2.5)
ax.grid(axis="x", color=GRID, lw=0.7, zorder=0); ax.set_axisbelow(True)
ax.set_title("The effect decays as the method tightens", fontsize=9, color=MUTED, pad=9)
save(fig, "fig4_hbond_decay")

# ---------------------------------------------------------------- F5 adducts
names = ["Native\nArg143", "MG-H1", "Argpyrimidine", "CEA"]
mean = [2.00, 1.52, 0.87, 2.43]
sd = [0.34, 0.14, 0.31, 0.38]
cols = [NATIVE, GLYC, GLYC, "#009E73"]
fig, ax = plt.subplots(figsize=(4.6, 3.2))
b = ax.bar(names, mean, yerr=sd, capsize=4, width=0.62, color=cols,
           edgecolor="white", linewidth=1.4, zorder=3,
           error_kw=dict(ecolor=MUTED, lw=1.2))
ax.axhline(mean[0], color=MUTED, lw=1, ls=(0, (4, 3)), zorder=2)
for i, (m, s_) in enumerate(zip(mean, sd)):
    if i:
        ax.text(i, m + s_ + 0.13, f"{m-mean[0]:+.2f}", ha="center", fontsize=9,
                color=INK, weight="bold")
ax.set_ylabel("Site-143 protein H-bonds")
ax.set_ylim(0, 3.15)
ax.grid(axis="y", color=GRID, lw=0.7, zorder=0); ax.set_axisbelow(True)
ax.set_title("CEA reverses direction — guanidinium retained", fontsize=9, color=MUTED, pad=9)
save(fig, "fig5_adduct_series")

print("\nAll five written to " + os.path.abspath(OUT))
