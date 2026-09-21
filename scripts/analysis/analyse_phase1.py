#!/usr/bin/env python3
"""
Phase 1: native vs MG-H1 monomer, explicit OPC solvent, n = 10, 300 ns.

STATISTICAL UNIT IS THE REPLICATE, never the frame. Frames within a trajectory
are autocorrelated; treating them as independent inflates n by ~10^4 and
manufactures significance. Each replicate contributes ONE number per endpoint.

With n = 10 per group the exact permutation test enumerates all C(20,10) =
184,756 splits, giving a floor of 1.1e-05. At the n = 5 of the first Phase 1
batch the floor was 0.0079; at the n = 3 of the implicit-solvent work it was
0.10, so nothing there could have reached significance.

Matched window: the final 150 ns of every replicate. Comparing different
windows across replicates is how this project has previously produced wrong
answers, so the window is fixed here and asserted.
"""
import os, sys, itertools, math
import numpy as np
import mdtraj as md

D = os.environ.get("MGO_DIR", "/scratch/a2148a01")
STRIDE = 1               # full.dcd is already 200 ps/frame
WINDOW_NS = 150.0        # final 150 ns of each 300 ns run
FRAME_PS = 200.0 * STRIDE
SITE = 142               # Arg143 / MGH143, zero-based
GLY61, CYS57 = 60, 56
HIS48 = 47               # third member of the CORE composite
CONTACT_CUT = 0.45        # nm; 4.5 A contact criterion
ALPHA_CORRECTED = 0.0042    # Bonferroni, 0.05 / 12 endpoints
LO, HI = 0.42, 0.50      # contact hysteresis band, nm (as in the implicit work)
N_BOOT = 10000
RNG = np.random.default_rng(20260904)

GROUPS = {"native": [f"native_mono_rep{i}" for i in range(1, 11)],
          "glyc":   [f"glyc_mono_rep{i}"   for i in range(1, 11)]}


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


def _resolve(tag, base):
    """Topology and trajectory path for TAG, in either layout.

    working:  $MGO_DIR/<base>.prmtop        + $MGO_DIR/prod/<tag>/full.dcd
    deposit:  $MGO_DIR/<base>_reduced.pdb   + $MGO_DIR/<tag>.dcd

    full.dcd holds every atom including water; the deposited trajectory holds the
    protein only, so the solvated prmtop does not match it and the reduced pdb
    written beside it does. Both carry MGH at residue 142, which prot.dcd did not.
    """
    pt = os.path.join(D, base + ".prmtop")
    pd = os.path.join(D, "prod", tag, "full.dcd")
    if os.path.exists(pt) and os.path.exists(pd):
        return md.load_prmtop(pt), pd
    rp = os.path.join(D, base + "_reduced.pdb")
    rd = os.path.join(D, tag + ".dcd")
    if os.path.exists(rp) and os.path.exists(rd):
        return md.load(rp).topology, rd
    return None, None


def analyse(tag):
    """Load from full.dcd, NOT prot.dcd.

    prot.dcd was written with select("protein"), and MDTraj does not classify
    MGH as protein -- so the adduct was absent from the glycated trajectories
    and every residue after it shifted by one. Index 142 pointed at LEU143 in
    the glycated systems and ARG142 in the native ones, i.e. two different
    residues were compared. full.dcd holds every atom, so it is immune.
    """
    base = tag.rsplit("_rep", 1)[0]
    top, dcd = _resolve(tag, base)
    if top is None:
        print("  !! no trajectory for " + tag + " under " + D); return None

    # fail loudly rather than silently measure the wrong residue again
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
        core_his=float((d_his < CONTACT_CUT).mean() * 100),
        core_cys=float((d_cys < CONTACT_CUT).mean() * 100),
        core_gly=float((d_gly < CONTACT_CUT).mean() * 100),
        core_occ=float(np.mean([(d_his < CONTACT_CUT).mean(),
                                (d_cys < CONTACT_CUT).mean(),
                                (d_gly < CONTACT_CUT).mean()]) * 100),
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
      f"window = final {WINDOW_NS:.0f} ns | permutation floor p = "
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
print("NOTE: switch rates come from 200 ps frames, so they are comparable")
print("BETWEEN these groups but NOT with the implicit-solvent switch rates,")
print("which were sampled far more finely.")
print("=" * 104)


print()
print("=" * 104)
print("PRIMARY ENDPOINT")
print("Site-143 CORE contact-network integrity: mean occupancy over His48/Cys57/Gly61,")
print(f"heavy-atom contact < {CONTACT_CUT*10:.1f} A.  Bonferroni alpha = {ALPHA_CORRECTED}.")
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
        _v = "  <-- SIGNIFICANT" if _p <= ALPHA_CORRECTED else "  <-- not significant"
    print(f"{_name:30s} {_a.mean():8.3f}+-{_a.std(ddof=1):5.3f} "
          f"{_b.mean():8.3f}+-{_b.std(ddof=1):5.3f} {_b.mean()-_a.mean():+9.3f} "
          f"[{_lo:+7.3f},{_hi:+7.3f}] {_p:8.4f} {_sep:>4s}{_v}")
print()
print("RANGE SEPARATION is reported in the 'sep' column.")
print("Every endpoint in the table above this block is EXPLORATORY,")
print("including Gly61 occupancy in isolation, which enters the analysis")
print("only as one component of the CORE composite.")
print("=" * 104)
