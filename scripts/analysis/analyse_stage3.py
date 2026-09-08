#!/usr/bin/env python3
"""
Stage 3: native vs MG-H1 apo-SH monomer, IMPLICIT GB-OBC2 solvent, n = 3, 300 ns.

This is the implicit-solvent arm of the solvent-model comparison. It is a port of
analyse_phase1.py (explicit OPC, n = 10) and deliberately runs the SAME endpoints,
the SAME window, and the SAME statistics, so that any difference between the two
arms is attributable to the solvent model rather than to the analysis.

STATISTICAL UNIT IS THE REPLICATE, never the frame. Frames within a trajectory
are autocorrelated; treating them as independent inflates n by ~10^4 and
manufactures significance. Each replicate contributes ONE number per endpoint.

With n = 3 per group the exact permutation test enumerates all C(6,3) = 20
splits, giving a floor of 0.10. NO ENDPOINT IN THIS ARM CAN REACH SIGNIFICANCE
AT ANY CONVENTIONAL ALPHA. That is a property of the design, not of the data:
read the effect sizes and confidence intervals, and ignore the p column except
as a reminder of the floor. The explicit arm at n = 10 has a floor of 1.1e-05.

FRAME RESOLUTION. These trajectories were written at 10 ps/frame; the explicit
runs were written at 200 ps/frame. Occupancies and H-bond counts are time
averages and are insensitive to this, but SWITCH RATES ARE NOT -- a finely
sampled trajectory reports more transitions for identical physics. The
trajectories are therefore strided by 20 to a matched 200 ps/frame, so all eight
endpoints are directly comparable between arms rather than only six of them.

Matched window: the final 150 ns of every replicate. Comparing different
windows across replicates is how this project has previously produced wrong
answers, so the window is fixed here and asserted.
"""
import os, sys, itertools, math
import numpy as np
import mdtraj as md

D = os.environ.get("MGO_DIR", "/home/jyoti/Projects/MGO")
STRIDE = 20              # traj.dcd is 10 ps/frame -> 200 ps, matching the explicit arm
WINDOW_NS = 150.0        # final 150 ns of each 300 ns run
FRAME_PS = 10.0 * STRIDE
SITE = 142               # Arg143 / MGH143, zero-based
GLY61, CYS57 = 60, 56
HIS48 = 47               # third member of the preregistered CORE
PREREG_CUT = 0.45        # nm; 4.5 A, the preregistered contact criterion
PREREG_ALPHA = 0.0042    # Bonferroni, fixed in PREREGISTRATION.md 2026-08-22
LO, HI = 0.42, 0.50      # contact hysteresis band, nm (as in the implicit work)
N_BOOT = 10000
RNG = np.random.default_rng(20260904)   # same seed as the explicit arm

GROUPS = {"native": [f"nativeSH_rep{i}" for i in range(1, 4)],
          "glyc":   [f"glycSH_rep{i}"   for i in range(1, 4)]}

# tag prefix -> topology, because the local runs do not follow the
# "<base>_rep<N>" / "<base>.prmtop" convention the explicit arm used.
TOPS = {"nativeSH": os.path.join(D, "sod1_new", "native_SH.prmtop"),
        "glycSH":   os.path.join(D, "sod1_new", "glyc_SH.prmtop")}


def switch_rate(d_nm, frame_ps):
    """Contacts formed+broken per 100 ns, with a hysteresis band so that
    thermal jitter across a single cutoff is not counted as a transition."""
    state, n = d_nm[0] < LO, 0
    for x in d_nm[1:]:
        if state and x > HI:
            state, n = False, n + 1
        elif not state and x < LO:
            state, n = True, n + 1
    return n / (len(d_nm) * frame_ps / 1e5)


def analyse(tag):
    """Implicit-solvent runs have no water, so traj.dcd holds the whole system
    and the prot.dcd/full.dcd hazard of the explicit arm does not arise here.
    The residue assertions are kept regardless -- they are what caught that bug."""
    base = tag.rsplit("_rep", 1)[0]
    top = md.load_prmtop(TOPS[base])
    dcd = os.path.join(D, tag, "traj.dcd")
    if not os.path.exists(dcd):
        print("  !! missing " + dcd); return None

    # fail loudly rather than silently measure the wrong residue
    assert top.residue(SITE).name in ("ARG", "MGH"), tag + " site is " + top.residue(SITE).name
    assert top.residue(GLY61).name == "GLY", tag + " GLY61 is " + top.residue(GLY61).name
    assert top.residue(CYS57).name == "CYS", tag + " CYS57 is " + top.residue(CYS57).name
    assert top.residue(HIS48).name in ("HIS", "HIE", "HID", "HIP"), \
        tag + " HIS48 is " + top.residue(HIS48).name

    t = md.load(dcd, top=top, stride=STRIDE)
    keep = int(WINDOW_NS * 1000 / FRAME_PS)
    assert t.n_frames >= keep, tag + ": too few frames"
    t = t[-keep:]

    d_gly = md.compute_contacts(t, [[SITE, GLY61]], scheme="closest-heavy")[0][:, 0]
    d_cys = md.compute_contacts(t, [[SITE, CYS57]], scheme="closest-heavy")[0][:, 0]
    d_his = md.compute_contacts(t, [[SITE, HIS48]], scheme="closest-heavy")[0][:, 0]

    at = t.topology.atom
    sub = t[::5]
    per_frame = [sum(1 for don, _, acc in fr
                     if at(don).residue.index == SITE
                     or at(acc).residue.index == SITE)
                 for fr in md.wernet_nilsson(sub, periodic=False)]

    # Global metrics over a COMMON residue set: all protein residues except
    # the site. Native has 153 protein residues, glycated 152 (MGH excluded),
    # so dropping the site leaves 152 in both and the numbers are strictly
    # comparable rather than differing by one residue.
    common = t.topology.select("protein and not resid " + str(SITE))
    prot = t.atom_slice(common)
    bb = prot.topology.select("backbone")
    ss = md.compute_dssp(prot, simplified=True)

    return dict(
        hbond_site=float(np.mean(per_frame)),
        gly61_occ=float((d_gly < LO).mean() * 100),
        gly61_switch=float(switch_rate(d_gly, FRAME_PS)),
        cys57_occ=float((d_cys < LO).mean() * 100),
        cys57_switch=float(switch_rate(d_cys, FRAME_PS)),
        rmsd=float(md.rmsd(prot, prot, 0, atom_indices=bb).mean() * 10),
        rg=float(md.compute_rg(prot).mean() * 10),
        beta=float((ss == "E").mean() * 100),
        core_his=float((d_his < PREREG_CUT).mean() * 100),
        core_cys=float((d_cys < PREREG_CUT).mean() * 100),
        core_gly=float((d_gly < PREREG_CUT).mean() * 100),
        core_occ=float(np.mean([(d_his < PREREG_CUT).mean(),
                                (d_cys < PREREG_CUT).mean(),
                                (d_gly < PREREG_CUT).mean()]) * 100),
        n_frames=t.n_frames, site=top.residue(SITE).name)


def perm_p(a, b):
    """Exact two-sided permutation test on the difference of means."""
    obs = abs(np.mean(b) - np.mean(a))
    pool = np.concatenate([a, b])
    k, hits, tot = len(a), 0, 0
    for idx in itertools.combinations(range(len(pool)), k):
        m = np.zeros(len(pool), bool); m[list(idx)] = True
        if abs(pool[~m].mean() - pool[m].mean()) >= obs - 1e-12:
            hits += 1
        tot += 1
    return hits / tot


def boot_ci(a, b):
    d = [np.mean(RNG.choice(b, len(b), True)) - np.mean(RNG.choice(a, len(a), True))
         for _ in range(N_BOOT)]
    return np.percentile(d, 2.5), np.percentile(d, 97.5)


print(__doc__)
res = {g: [] for g in GROUPS}
for g, tags in GROUPS.items():
    for tag in tags:
        r = analyse(tag)
        if r:
            res[g].append(r)
            print(f"  {tag:20s} frames={r['n_frames']:6d}  "
                  f"site={r['site']}  hb={r['hbond_site']:5.2f}  "
                  f"gly61={r['gly61_occ']:5.1f}%  "
                  f"cys57={r['cys57_occ']:5.1f}%", flush=True)

METRICS = [("hbond_site",   "Site-143 H-bonds (protein)"),
           ("gly61_occ",    "Site-Gly61 occupancy (%)"),
           ("gly61_switch", "Site-Gly61 switches /100 ns"),
           ("cys57_occ",    "Site-Cys57 occupancy (%)"),
           ("cys57_switch", "Site-Cys57 switches /100 ns"),
           ("rmsd",         "Backbone RMSD (A)"),
           ("rg",           "Radius of gyration (A)"),
           ("beta",         "Beta-sheet content (%)")]

print("\n" + "=" * 104)
print(f"n = {len(res['native'])} native, {len(res['glyc'])} glycated | "
      f"window = final {WINDOW_NS:.0f} ns | frame = {FRAME_PS:.0f} ps | "
      f"permutation floor p = "
      f"{2/math.comb(len(res['native'])+len(res['glyc']), len(res['native'])):.3g}")
print("=" * 104)
print(f"{'endpoint':30s} {'native':>16s} {'glycated':>16s} {'delta':>9s} "
      f"{'95% CI':>20s} {'p':>8s} {'sep':>4s}")
rows = []
for key, name in METRICS:
    a = np.array([r[key] for r in res["native"]])
    b = np.array([r[key] for r in res["glyc"]])
    if len(a) < 2 or len(b) < 2:
        continue
    lo, hi = boot_ci(a, b)
    p = perm_p(a, b)
    sep = "yes" if (max(a) < min(b) or max(b) < min(a)) else "no"
    rows.append((name, p))
    print(f"{name:30s} {a.mean():8.3f}+-{a.std(ddof=1):5.3f} "
          f"{b.mean():8.3f}+-{b.std(ddof=1):5.3f} {b.mean()-a.mean():+9.3f} "
          f"[{lo:+7.3f},{hi:+7.3f}] {p:8.4f} {sep:>4s}")

# Benjamini-Hochberg, STEP-UP: find the largest passing rank, then every
# rank at or below it is significant. Marking rows independently wrongly
# failed rank 1 while passing rank 2 at an identical p-value.
m = len(rows)
order = sorted(range(m), key=lambda i: rows[i][1])
kmax = 0
for rank, i in enumerate(order, 1):
    if rows[i][1] <= 0.05 * rank / m:
        kmax = rank
print()
print("Benjamini-Hochberg (FDR 0.05), step-up:")
for rank, i in enumerate(order, 1):
    thr = 0.05 * rank / m
    mark = "SIGNIFICANT" if rank <= kmax else ""
    print(f"  {rank:2d}. {rows[i][0]:32s} p={rows[i][1]:.4f}  thr={thr:.4f}  {mark}")
print()
print(f"  {kmax} endpoint(s) significant at FDR 0.05" if kmax
      else "  no endpoint significant at FDR 0.05")
print()
print("NOTE: n = 3 per group puts the permutation floor at 0.10, so the p column")
print("here cannot fall below that no matter how large the effect. Compare the")
print("DELTAS and CONFIDENCE INTERVALS against the explicit arm, not the p-values.")
print("Switch rates ARE comparable between arms: both are computed at 200 ps.")
print("=" * 104)


# ---------------------------------------------------------------------------
# PREREGISTERED PRIMARY ENDPOINT
# Reported separately from the exploratory BH family above, and against the
# Bonferroni alpha fixed in PREREGISTRATION.md, because that is what was
# promised. Folding it into the BH family would move the thresholds for the
# eight endpoints already reported.
# ---------------------------------------------------------------------------
print()
print("=" * 104)
print("PREREGISTERED PRIMARY ENDPOINT  (PREREGISTRATION.md, sealed 2026-08-22 00:30 KST)")
print(f"Site-143 CORE contact-network integrity: mean occupancy over His48/Cys57/Gly61,")
print(f"heavy-atom contact < {PREREG_CUT*10:.1f} A.  Preregistered Bonferroni alpha = {PREREG_ALPHA}.")
print("=" * 104)
print(f"{'endpoint':30s} {'native':>16s} {'glycated':>16s} {'delta':>9s} "
      f"{'95% CI':>20s} {'p':>8s} {'sep':>4s}")
for _k, _name in [("core_occ", "CORE composite (PRIMARY)"),
                  ("core_his", "   component His48"),
                  ("core_cys", "   component Cys57"),
                  ("core_gly", "   component Gly61")]:
    _a = np.array([r[_k] for r in res["native"]])
    _b = np.array([r[_k] for r in res["glyc"]])
    if len(_a) < 2 or len(_b) < 2:
        continue
    _lo, _hi = boot_ci(_a, _b)
    _p = perm_p(_a, _b)
    _sep = "yes" if (max(_a) < min(_b) or max(_b) < min(_a)) else "no"
    _v = ""
    if _k == "core_occ":
        _v = "  <-- SIGNIFICANT (prereg)" if _p <= PREREG_ALPHA else "  <-- not significant (prereg)"
    print(f"{_name:30s} {_a.mean():8.3f}+-{_a.std(ddof=1):5.3f} "
          f"{_b.mean():8.3f}+-{_b.std(ddof=1):5.3f} {_b.mean()-_a.mean():+9.3f} "
          f"[{_lo:+7.3f},{_hi:+7.3f}] {_p:8.4f} {_sep:>4s}{_v}")
print()
print("RANGE SEPARATION was the preregistered PRIMARY EVIDENCE ('sep' column);")
print("permutation p and the bootstrap CI are reported alongside it. Note that")
print("the plan specified Welch p as descriptive -- the exact permutation test")
print("used here is a DEVIATION, and a more conservative one at this n.")
print()
print("Every endpoint in the table above this block is EXPLORATORY / POST HOC")
print("with respect to PREREGISTRATION.md, including Gly61 occupancy in")
print("isolation, which the plan named only as one component of the composite.")
print("=" * 104)
