#!/usr/bin/env python3
"""P6 -- no sequential-damage propagation (NULL prediction).

Prediction: modifying Arg143 will produce REDISTRIBUTION of solvent accessibility
among the other arginine and lysine side chains, not a systematic increase.
Refuted if distal Arg/Lys exposure increases systematically and survives
multiple-comparison correction.

Basis (from PREREGISTRATION.md): on the completed apo-SS data three distal
lysines exceeded replicate scatter but moved in OPPOSITE directions (Lys122
-29%, Lys136 +26%, Lys128 +13%), and nothing survived correction across 14
sites. This restates that as a prediction for the apo-SH condition, which now
has matched 300 ns data.

Design. Per-residue solvent-accessible surface area (Shrake-Rupley) over the
final 150 ns, native vs glycated, n = 3 per group. The panel is every ARG and
LYS in the native topology EXCEPT the modified site, so the two conditions are
compared over an identical residue set. Exact permutation per residue,
Benjamini-Hochberg across the panel at FDR 0.05.

n = 3 per group: the permutation floor is 0.10, so no residue can reach
significance. The prediction concerns DIRECTION -- whether the changes share a
sign -- which is testable at this n by a sign test across the panel even though
the individual comparisons are not.
"""
import os, sys, itertools, math
import numpy as np
import mdtraj as md

D = os.environ.get("MGO_DIR", "/home/jyoti/Projects/MGO")
STRIDE = 20
FRAME_PS = 10.0 * STRIDE
WINDOW_NS = 150.0
SITE = 142
SASA_STRIDE = 5          # SASA is expensive; 150 frames over the window
N_BOOT = 10000
RNG = np.random.default_rng(20260904)

TOP_N = f"{D}/sod1_new/native_SH.prmtop"
TOP_G = f"{D}/sod1_new/glyc_SH.prmtop"
GROUPS = {"native": (TOP_N, [f"nativeSH_rep{i}" for i in (1, 2, 3)]),
          "glyc":   (TOP_G, [f"glycSH_rep{i}"   for i in (1, 2, 3)])}


def panel(top):
    """Every ARG/LYS except the modified site, by residue index."""
    return [r.index for r in top.residues
            if r.name in ("ARG", "LYS") and r.index != SITE]


def sasa_profile(top_f, tag, idx):
    top = md.load_prmtop(top_f)
    dcd = f"{D}/{tag}/traj.dcd"
    if not os.path.exists(dcd):
        print(f"  !! missing {dcd}"); return None
    t = md.load(dcd, top=top, stride=STRIDE)
    keep = int(WINDOW_NS * 1000 / FRAME_PS)
    assert t.n_frames >= keep, f"{tag}: {t.n_frames} frames, need {keep}"
    t = t[-keep:][::SASA_STRIDE]
    sasa = md.shrake_rupley(t, mode="residue")      # nm^2, (frames, residues)
    return sasa[:, idx].mean(axis=0) * 100.0        # -> A^2


def perm_p(a, b):
    obs = abs(np.mean(b) - np.mean(a))
    pool = np.concatenate([a, b]); k = len(a); hits = tot = 0
    for i in itertools.combinations(range(len(pool)), k):
        m = np.zeros(len(pool), bool); m[list(i)] = True
        if abs(pool[~m].mean() - pool[m].mean()) >= obs - 1e-12:
            hits += 1
        tot += 1
    return hits / tot


def boot_ci(a, b):
    d = [np.mean(RNG.choice(b, len(b), True)) - np.mean(RNG.choice(a, len(a), True))
         for _ in range(N_BOOT)]
    return np.percentile(d, 2.5), np.percentile(d, 97.5)


print(__doc__)
top_n = md.load_prmtop(TOP_N)
idx = panel(top_n)
labels = [f"{top_n.residue(i).name}{top_n.residue(i).resSeq + 1}" for i in idx]
print(f"panel: {len(idx)} residues (site index {SITE} excluded)")
print("  " + "  ".join(labels))

prof = {}
for g, (top_f, tags) in GROUPS.items():
    prof[g] = []
    for tag in tags:
        v = sasa_profile(top_f, tag, idx)
        if v is not None:
            prof[g].append(v)
            print(f"  {tag:16s} mean panel SASA = {v.mean():6.2f} A^2", flush=True)

A = np.array(prof["native"])       # (3, n_res)
B = np.array(prof["glyc"])
floor = 2 / math.comb(len(A) + len(B), len(A))

print("\n" + "=" * 92)
print(f"P6 -- distal Arg/Lys SASA, apo-SH, final {WINDOW_NS:.0f} ns, "
      f"n = {len(A)} per group, floor p = {floor:.2f}")
print("=" * 92)
print(f"{'residue':>10s} {'native':>14s} {'glycated':>14s} {'delta':>8s} "
      f"{'%':>8s} {'95% CI':>18s} {'p':>7s}")
rows = []
for j, lab in enumerate(labels):
    a, b = A[:, j], B[:, j]
    d = b.mean() - a.mean()
    pct = 100 * d / a.mean() if a.mean() else float("nan")
    lo, hi = boot_ci(a, b); p = perm_p(a, b)
    rows.append((lab, d, pct, p))
    print(f"{lab:>10s} {a.mean():8.2f}+-{a.std(ddof=1):4.2f} "
          f"{b.mean():8.2f}+-{b.std(ddof=1):4.2f} {d:+8.2f} {pct:+7.1f}% "
          f"[{lo:+7.2f},{hi:+7.2f}] {p:7.2f}")

# Benjamini-Hochberg, step-up
m = len(rows)
order = sorted(range(m), key=lambda i: rows[i][3])
kmax = 0
for rank, i in enumerate(order, 1):
    if rows[i][3] <= 0.05 * rank / m:
        kmax = rank
print(f"\nBenjamini-Hochberg (FDR 0.05) across {m} sites: "
      f"{kmax} significant" if kmax else
      f"\nBenjamini-Hochberg (FDR 0.05) across {m} sites: none significant")

# direction test -- this is what P6 is actually about
up = sum(1 for _, d, _, _ in rows if d > 0)
dn = m - up
k = max(up, dn)
sign_p = 2 * sum(math.comb(m, i) for i in range(k, m + 1)) / 2 ** m
sign_p = min(1.0, sign_p)
print("\n" + "=" * 92)
print(f"DIRECTION: {up} of {m} sites increase, {dn} decrease.  "
      f"Two-sided sign test p = {sign_p:.3f}")
big = sorted(rows, key=lambda r: -abs(r[2]))[:4]
print("largest fractional changes: " +
      ", ".join(f"{lab} {pct:+.0f}%" for lab, _, pct, _ in big))
if sign_p > 0.05:
    print("\n--> Changes do NOT share a systematic direction, and no site survives")
    print("    correction. P6 (redistribution, not systematic increase) SUPPORTED.")
else:
    print("\n--> Changes DO share a direction. Check whether any site survives")
    print("    correction before calling P6 refuted.")
print("=" * 92)
